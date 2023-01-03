#!/bin/bash

data="BRCA"
sh_path=$(pwd)
echo "  "
echo "-----------------------------------------------------------------------"
echo "                                       cancer type: $data                                       "
echo "-----------------------------------------------------------------------"
echo "  "
echo "Creating filefolders for each output..."
echo "-----------------------------------------------------------------------"
echo "  "
mkdir -p ../PPR_data/"$data"/
mkdir -p ./getpairs_sh/
for str1 in "hi" "iref" "mult"
do
    mkdir -p ../sub_mat/"$data"/"$str1"/
    for i in 8
    do
        for j in 0 1
        do
            for str2 in "0.1"
            do
                rm -rf ../out_module/"$data"/"$str1"/e"$str2"_l"$j"_u"$i"
                mkdir -p ../out_module/"$data"/"$str1"/e"$str2"_l"$j"_u"$i"/
                rm -rf ../CLK_module/"$data"/"$str1"/e"$str2"_l"$j"_u"$i"
                mkdir -p ../CLK_module/"$data"/"$str1"/e"$str2"_l"$j"_u"$i"/
            done
        done
    done
done
for str3 in "0.1"
do
    mkdir -p ../genepairs/e"$str3"_"$data"/
    mkdir -p ../consensus/e"$str3"_"$data"/
done

echo "  "
echo "Processing data..."
echo "-----------------------------------------------------------------------"
echo "  "
for str1 in "hi" "iref" "mult"
do
{
    echo "$str1"
    echo "  "
    time   ../src/prep_data ../mutation_data/"$data" ../PPR_data/PPR_"$str1".m ../PPR_data/"$str1"_index_file.txt 15 ../PPR_data/"$data"/"$str1"_
    echo "  "
} &
done
wait

echo "  "
echo "Creating the mutual matrix for each modules..."
echo "-----------------------------------------------------------------------"
echo "  "
for str1 in "hi" "iref" "mult" 
do
{
    echo "$str1"
    echo "  "
    time   ../src/creat_subnet ../PPR_data/"$data"/"$str1"_C ../PPR_data/"$data"/"$str1"_Matrix01 300 2 5 ../sub_mat/"$data"/"$str1"/
    echo "  "
} &
done
wait

echo "  "
echo "Selecting random networks..."
echo "-----------------------------------------------------------------------"
echo "  "
for str1 in "hi" "iref" "mult"
do
{
    echo "$str1"
    echo "  "
    time   ../src/entropy_01 ../PPR_data/"$data"/"$str1"_C ../PPR_data/"$data"/"$str1"_Matrix01 2 5 100000 1 ../PPR_data/"$data"/"$str1"_
echo "  "
} &
done
wait

echo "  "
echo "Running MaxCLK programing..."
echo "-----------------------------------------------------------------------"
echo "  "
#time
thread_num=6
tmpfile="temp_fifo"
mkfifo -m 777 $tmpfile
exec 6<>${tmpfile}
for i in 1 2 3 4 5 6
do
{
    echo "" >&6
}
done

for str1 in hi iref mult
do
    Gene_num=$(cat ../PPR_data/"$data"/"$str1"_Gene_num.txt)
    echo "$data $str1 Number of mutation genes: $Gene_num ..."
    echo "  "
    for i in 8
    do
        for j in 0 1
        do
            for nset in 2 3 4 5
            do
                for str2 in "0.1"
                do
                    echo "Run MaxCLK with ... $data $str1 __ Size of clique: $nset errEx: $j MMN: $i ..."
                    echo " "
                    for s in `seq 1 1 $Gene_num`
		    do
                    {
                        read -u 6
                        {
                            if [ -s ../sub_mat/"$data"/"$str1"/mMat"$s"_t"$nset".txt ]; then
                            echo $data $str1 $str2 $i $j $s $nset >> ../null
                            ../src/maxCliques ../sub_mat/"$data"/"$str1"/mMat"$s"_t"$nset".txt $nset $i $j $str2 ../CLK_module/"$data"/"$str1"/e"$str2"_l"$j"_u"$i"/CLKmodule_"$nset".txt >> ../null
			    fi
                            echo "" >&6
		        } &
                    }
                    done
                    wait
                    ../src/filterCLK ../CLK_module/"$data"/"$str1"/e"$str2"_l"$j"_u"$i"/CLKmodule_"$nset".txt ../PPR_data/"$data"/"$str1"_Considered_genes ../PPR_data/"$data"/"$str1"_Ex_entropy.txt 0.8 0.05 ../out_module/"$data"/"$str1"/e"$str2"_l"$j"_u"$i"/SKmodule_"$nset".txt
                    rm ../CLK_module/"$data"/"$str1"/e"$str2"_l"$j"_u"$i"/CLKmodule_"$nset".txt
                    ../src/minset_greedy ../PPR_data/"$data"/"$str1"_Matrix01 ../out_module/"$data"/"$str1"/e"$str2"_l"$j"_u"$i"/SKmodule_"$nset".txt ../PPR_data/"$data"/"$str1"_Considered_genes ../out_module/"$data"/"$str1"/e"$str2"_l"$j"_u"$i"/CSKmodule_"$nset".txt ../out_module/"$data"/"$str1"/e"$str2"_l"$j"_u"$i"/CSKmodule_"$nset"_genes.txt
                    rm ../out_module/"$data"/"$str1"/e"$str2"_l"$j"_u"$i"/SKmodule_"$nset".txt
                done
            done
        done
    done
rm -rf $sh_path/../sub_mat/"$data"/"$str1"/
done
exec 6>&-
rm ${tmpfile}

echo "  "
echo "Running consensus programing..."
echo "-----------------------------------------------------------------------"
echo "  "
#time   
       for s in 1 2 3 4 5 6 7
       do
           for str2 in "0.1"
           do
               echo "#!/bin/bash/" > ./getpairs_sh/getpairs_"$data"_"$str2"_"$s".sh
               echo "$sh_path/../src/getgenepairs $sh_path/../genepairs/e"$str2"_"$data"/"$s"_ "$s" \\" >> ./getpairs_sh/getpairs_"$data"_"$str2"_"$s".sh
               for j in 0 1
               do
                   for i in 8
                   do
                       for m in 2 3 4 5
                       do
                           for str1 in "hi" "iref" "mult"
                           do
                               echo " $sh_path/../out_module/"$data"/"$str1"/e"$str2"_l"$j"_u"$i"/CSKmodule_"$m".txt \\" >> ./getpairs_sh/getpairs_"$data"_"$str2"_"$s".sh
                           done
                       done
                   done
               done
           done
       done
#time   
       for str2 in "0.1"
       do
           for s in 1 2 3 4 5 6 7
           do
               sh ./getpairs_sh/getpairs_"$data"_"$str2"_"$s".sh
           done
       done
#time
       for str2 in "0.1"
       do
           for s in 1 2 3 4 5 6 7
           do
               echo "$data e=$str2 s=$s"
               ../src/consensus ../genepairs/e"$str2"_"$data"/"$s"_SWgenepairs.txt ../consensus/e"$str2"_"$data"/subnet_A_"$s".txt ../consensus/e"$str2"_"$data"/subnet_B_"$s".txt 30 ../ref_genes/NCGgenes.txt
           done
       done
echo "  "
echo "Tips:   "
echo "  "
echo "-----------------------------------------------------------------------"
echo "Congratulations!"
echo "-----------------------------------------------------------------------"
echo "  "

