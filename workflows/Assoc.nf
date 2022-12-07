//import processes
include { option_check } from '../modules/option_check.nf'
include { prefilter } from '../modules/prefilter.nf'
include { extract_dosage } from '../modules/extract_dosage.nf'
include { generate_pcs } from '../modules/flashpca.nf'
include { make_saige_covars } from '../modules/make_saige_covars.nf'

include { plink_assoc;
          merge_plink_lifted
        } from '../modules/plink1.nf'

include { prune;
		  merge_plink;
          make_plink } from '../modules/plink2.nf'

include { liftover_pruned;
          liftover 
		} from '../modules/liftover.nf'

include { gen_r2_list;
          split_vcf_ranges } from '../modules/bcftools.nf'

include { saige_null; 
          saige_assoc
        } from '../modules/saige.nf'

include { merge_saige_results; 
          merge_liftover_tables;
          merge_plink_results;
          merge_r2
        } from '../modules/merge_processes.nf'
		
include {regenie_step1;
		 regenie_step2;
		 phenofile_from_fam } from '../modules/regenie.nf'

include { prune_python_helper } from '../modules/python2.nf'
include { lift_plink_sumstats;
		//lift_regenie_sumstats;
		  lift_saige_sumstats } from '../modules/liftsumstats.nf'

//function definitions
def get_file_details(filename) {
    def m = filename =~ /\/([^\/]+)(\.vcf\.gz|\.bgen)$/
    return [ m[0][1], m[0][2] ]
	}

def get_chromosome_code(filename) {
    def m = filename =~ /\/([^\/]+).ap_prf.vcf.gz$/
    return m[0][1]
	}


def check_fam_for_saige(file) {
	def lines = file.readLines()
        int x = 0
    	int y = 0	
	lines.each { String line ->
  	if(line.split(" |\t")[5] == "2") x++
  	if(line.split(" |\t")[5] == "1") y++

	}
 	return(x > 50 && y > 50)

	}

//PIPELINE WORKFLOW
workflow assoc{

	main:

	option_check()
	
	//if(!params.fam){exit 1, "Cannot find fam file"}

	params.fam_length=check_fam_for_saige( file(params.fam) )

	for_saige_imp = Channel.fromPath(params.input_imputed_glob, checkIfExists: true).map { it ->
    	def match = get_file_details(it)
    	[it, match[1]] }

	prefilter( for_saige_imp,
				Channel.fromPath(params.fam) )


	ch_mapped_prefilter = prefilter.out.map { it -> [it[0], it[1], get_chromosome_code(it[0]), it[2]] }//.into { for_split_vcf; for_extract_dosage; for_gen_r2; for_plink_imp }

	extract_dosage( ch_mapped_prefilter )

	gen_r2_list( ch_mapped_prefilter )

	merge_r2( gen_r2_list.out.collect() )

	make_plink ( ch_mapped_prefilter,
				 Channel.fromPath(params.fam) )

	merge_plink ( make_plink.out.collect() )

	//removed python part from the original prune process
	prune_python_helper(merge_plink.out,
						merge_r2.out )
	prune( merge_plink.out,
		   merge_r2.out,
		   prune_python_helper.out )

//FLASHPCA2
	generate_pcs( prune.out )

	make_saige_covars( generate_pcs.out,
					   Channel.fromPath(params.fam, checkIfExists: true ).ifEmpty { exit 1, "Cannot find fam file"} )
//PLINK
	plink_assoc( extract_dosage.out.combine(make_saige_covars.out.for_plink ) )

	merge_plink_results( plink_assoc.out.collect() )

//LIFTOVER
	if(params.ucsc_liftover != "" && !params.disable_liftover){
		liftover( make_plink.out )

		liftover_pruned( prune.out )

		merge_plink_lifted( liftover.out.mergelifted )

		merge_liftover_tables( liftover.out.for_merge_tables.collect() )

		lift_plink_sumstats( merge_plink_results.out,
							 merge_liftover_tables.out )
	}

//SAIGE
	if(params.fam_length && params.run_saige){
		saige_null( prune.out, 
					make_saige_covars.out.covars )

		split_vcf_ranges( ch_mapped_prefilter )

		saige_assoc( split_vcf_ranges.out.transpose().combine(saige_null.out.toList()).map{ it -> [ it[0], it[1], it[2], it[3], it[4], it[5], it[6][0], it [6][1] ] },
					 Channel.fromPath(params.fam) )

		merge_saige_results( saige_assoc.out.collect(),
							 make_saige_covars.out.covars )

		//lift only when liftover was activated
		if(params.ucsc_liftover != "" && !params.disable_liftover){
			lift_saige_sumstats( merge_saige_results.out,
								 merge_liftover_tables.out )
		}
	}

//REGENIE
	if(!params.disable_regenie){
		if(params.phenofile){
			//TODO: check if phenofile exist
			ch_pheno = Channel.fromPath(params.phenofile, checkIfExists: true ).ifEmpty { exit 1, "Cannot find phenofile"}
		}else{
			phenofile_from_fam( Channel.fromPath(params.fam, checkIfExists: true ).ifEmpty { exit 1, "Cannot find fam file"} )
			ch_pheno = phenofile_from_fam.out
		}
		//Regenie step1 should be run with less than 1mio SNPs, therefor we use the pruned plink-files
		regenie_step1( prune.out,
					   make_saige_covars.out.covars,
					   ch_pheno )

		regenie_step2( merge_plink.out,
					   make_saige_covars.out.covars,
					   regenie_step1.out,
					   ch_pheno )
		//TODO: lift_regenie_sumstats
		/*
		if(params.ucsc_liftover != "" && !params.disable_liftover){
			lift_regenie_sumstats( regenie_step2.out,
								   merge_liftover_tables.out )
		}
		*/
	}
}
