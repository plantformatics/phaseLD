for i in `seq 1 9`;
do
	grep chr0$i $1 > chr0$i.gen
	sliding_window_LD.pl chr0$i.gen > chr0$i.ld
	filter_ld_snps.pl chr0$i.ld chr0$i.gen > chr0$i.ld.gen
done

for i in `seq 10 12`;
do
	grep chr$i $1 > chr$i.gen
	sliding_window_LD.pl chr$i.gen > chr$i.ld
	filter_ld_snps.pl chr$i.ld chr$i.gen > chr$i.ld.gen
done

cat *.ld.gen > $1.ld_filt.gen
