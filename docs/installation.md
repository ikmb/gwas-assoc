# Installation

## At the IKMB

If you are at the IKMB, you will not have to do anything to make this run, it is all pre-configured for our compute system(s).
 
## Site-specific config file

This pipeline requires a site-specific configuration file to be able to talk to your local cluster or compute infrastructure. Nextflow supports a wide
range of such infrastructures, including Slurm, LSF and SGE - but also Kubernetes and AWS. For more information, see [here](https://www.nextflow.io/docs/latest/executor.html).

Please see conf/medcluster.config for an example of how to configure this pipeline for a Slurm queue.

All software is provided through Docker containers - this requires for your compute system to run either Docker or Singularity (more common on HPC systems). Details on how to specify singularity as your container engine are provided in the config file for our medcluster (medcluster.config).

With this information in place, you will next have to create an new site-specific profile for your local environment in `nextflow.config` using the following format:

```

profiles {
	
	your_profile {
		includeConfig 'conf/base.config'
		includeConfig 'conf/your_cluster.config'
	}
}

```

This would add a new profile, called `your_profile` which uses (and expects) conda to provide all software. 

`base.config` Basic settings about resource usage for the individual pipeline stages. 

`resources.config` Gives information about the files that are to be used during analysis for the individual human genome assemblies.

