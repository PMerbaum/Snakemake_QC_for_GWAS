mkdir -p logs/slurm

snakemake \
-j450 \
--latency-wait=300 \
--cluster="sbatch -c {threads} \
           --mem={resources.mem_mb}M \
           --time={resources.time} \
           --output=logs/slurm/%j.out \
           --error=logs/slurm/%j.out \
           --parsable" \
--cluster-status ./config/cluster_status.py