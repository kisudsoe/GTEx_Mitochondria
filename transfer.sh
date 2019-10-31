# 191031
#scp -i SKim-EC2.pem src/bedtools_closest_roadmap.r src/regulome.r src/lncrnasnp.r ec2-user@ec2-54-89-92-92.compute-1.amazonaws.com:~/src
#scp -i SKim-EC2.pem db/encode_tfbs_merge.bed db/roadmap_enh_merge.bed db/lncRNASNP2_snplist.txt.rds ec2-user@ec2-54-89-92-92.compute-1.amazonaws.com:~/db
scp -i SKim-EC2.pem annotation.sh ec2-user@ec2-54-89-92-92.compute-1.amazonaws.com:~/