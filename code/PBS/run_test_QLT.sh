#!/bin/bash
for k in 20 50 100
#for k in 20
do
	for data_fn in generate_data_QLT generate_data_high_dimensional
	#for data_fn in generate_data_QLT
	do
		for start in cold_start warm_start_covariates warm_start_likelihood
		#for start_type in cold_start
		do
			#for prior in log_prob1 BIC ZE 3 4 5 6 7
			for prior in log_prob1 BIC ZE 3 4 5 6 7
			do
				echo $k $data_fn $start $prior
				export k=$k
				export data_fn=$data_fn
				export start_type=$start_type
				export prior=$prior
				#echo qsub -P RDS-FSC-ZIPVB-RW -q small ~/test_QLT.pbs
				qsub -P RDS-FSC-ZIPVB-RW -q small-express -v k,data_fn,start_type,prior ~/test_QLT_backup.pbs
			done
		done
	done
done
