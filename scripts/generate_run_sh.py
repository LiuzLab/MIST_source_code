import os

data_names = os.listdir("../data/")
data_names = [dn for dn in data_names if os.path.isdir(f'../data/{dn}')]

with open("run_holdout_experiments.sh", "w") as f:
	for dn in data_names:
		dp = os.path.abspath(f'../data/{dn}')
		f.write(f"nohup python3 run_holdout_experiments.py {dp} > ./logs/{dn}.out &\n")