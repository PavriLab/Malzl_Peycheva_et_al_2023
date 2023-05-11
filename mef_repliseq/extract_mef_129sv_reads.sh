for fq in `ls fastqs/*`;
do
	python3 utils/pyharp.py -f ${fq} -bti1 /groups/pavri/bioinfo/daniel/analysisMEF/repliseq/btindex/129_sv/129_sv -bti2 /groups/pavri/bioinfo/daniel/analysisMEF/repliseq/btindex/cast_ei/cast_ei -t 16 -o processed
done
