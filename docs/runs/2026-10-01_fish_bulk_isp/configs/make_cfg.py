import yaml, sys, copy
run = sys.argv[1]
cfg = yaml.safe_load(open("config/config.yaml"))
base = f"runs/{run}"
cfg["paths"] = {k: f"{base}/{v.split('/')[-1]}" for k, v in cfg["paths"].items()}
b = cfg["bulk_isp"]
b.update(counts="/mnt/nas/kaisar/te-oa-full/counts-singlebatch/all",
         samples="runs/fish_samples.tsv", stratify_by=["tissue", "treatment"],
         extra_regulator_patterns=["^znf438$"])
yaml.safe_dump(cfg, open(f"runs/{run}.yaml", "w"), sort_keys=False)
s = copy.deepcopy(cfg); s["paths"] = {k: v.replace(run, run + "_no_outliers") for k, v in cfg["paths"].items()}
s["bulk_isp"]["exclude_where"] = {"column": "source_study_outlier", "value": "yes"}
yaml.safe_dump(s, open(f"runs/{run}_no_outliers.yaml", "w"), sort_keys=False)
