#!/bin/bash
ulimit -s unlimited
start_time=`date +%s`
count=1
result=log/`date +%Y%m%d_%H-%M-%S`.log
echo "---a-expan ---"
for T in 1 2 3 4 5 6 7 8
do
    for lambda in 5
    do
        # range_size <= label_size
        job_start=`date +%s`
        ./a-estr input/tsukuba_ output/tsukuba_${lambda}_${T}.bmp 16 $lambda $T >>  ${result}
        job_end=`date +%s`
        time=$((job_end - job_start));
        count=`expr $count + 1`
        echo "tsukuba T=${T}, lambda=${lambda} [${time}s]";
    done
done
#
# for T in 1 2 3 4 5 6 7 8
# do
#     for lambda in 1 2 3 4 5 6 7
#     do
#         # range_size <= label_size
#         job_start=`date +%s`
#         ./a-estr input/venus_ output/venus_${lambda}_${T}.bmp 8 $lambda $T >>  ${result}
#         job_end=`date +%s`
#         time=$((job_end - job_start));
#         count=`expr $count + 1`
#         echo "venus T=${T}, lambda=${lambda} [${time}s]";
#     done
# done
#
# for T in 1 2 3 4 5 6 7 8
# do
#     for lambda in 1 2 3 4 5 6 7
#     do
#         # range_size <= label_size
#         job_start=`date +%s`
#         ./a-estr input/teddy_ output/teddy_${lambda}_${T}.bmp 4 $lambda $T  >>  ${result}
#         job_end=`date +%s`
#         time=$((job_end - job_start));
#         count=`expr $count + 1`
#         echo "teddy T=${T}, lambda=${lambda} [${time}s]";
#     done
# done

end_time=`date +%s`
time=$((end_time -start_time));
rm temp.txt

echo "@trsk_1st 全ての処理が完了しましたっ(a-expan)! 総所要時間[${time}s]" | tw --pipe --user="trsk_1st"

git add ${result}
git commit -m "job_${result}"
git push origin master
echo "------"
