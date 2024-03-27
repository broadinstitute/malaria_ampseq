#AMP\_SEQ

Check point 1

Quick git push: git add . && git commit -m "Updated wdl" && git push origin main

Quick docker push: 
#Comment out the last two lines of the file
docker build --platform linux/amd64 -t jorgeamaya/mixed_reads_ampseq . && docker push jorgeamaya/mixed_reads_ampseq

## Development

Build docker: docker build --platform linux/amd64 -t jorgeamaya/mixed_reads_ampseq .

Run docker container with mounted volume: docker run -v ~/Desktop/Data_Repository_AmpSeq/SIMPLseq_CI/Plate_1:/Data --platform linux/amd64 -it jorgeamaya/mixed_reads_ampseq bash

Quick build and exec: docker build --platform linux/amd64 -t jorgeamaya/mixed_reads_ampseq . && docker run -v ~/Desktop/Data_Repository_AmpSeq/SIMPLseq_CI/Plate_1:/Data --platform linux/amd64 -it jorgeamaya/mixed_reads_ampseq bash

Copy files if necessary: container_name=$(docker ps -qf "name=$1") && docker cp /Users/jar4142/Desktop/Data_Repository_AmpSeq/SIMPLseq_CI/Plate_1/barcodes_matches.csv "$container_name":barcodes_matches.csv && container_name=$(docker ps -qf "name=$1") && docker cp /Users/jar4142/Desktop/Plate_1_Results/missing_files.tsv "$container_name":missing_files.tsv

Copy directory back to host machine from running docker container: container_name=$(docker ps -qf "name=$1") && docker cp "$container_name":/Results .

If at any point the previous command doesn't work
1. Obtain the name of the container with the function: docker ps -q or from the prompt that points to the running container. For example, in the prompt "(ampseq_env) root@9f93f54a5455:/#" the name of the container is 9f93f54a5455.
2. Replace <container_name> with the name of the container and run the following command: docker cp <container_name>:Results .

Clearing Docker containers to realese memory: 
```
docker stop $(docker ps -a -q) && docker rm $(docker ps -a -q)
docker container prune; docker image prune; docker volume prune
```

Rscript /render_report.R -d "/Results/Merge/" -o "/Report/" -p "barcodes_matches.csv" -m 1000 -c 0.5 -mf "Results/missing_files.tsv"
