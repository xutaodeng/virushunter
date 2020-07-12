import os

file=["Sample_167_1_K27Ac_Target_ChIPseeqer.txt" \
,"Sample_167_2_K27Ac_Target_ChIPseeqer.txt" \
,"Sample_167_3_K27Ac_Target_ChIPseeqer.txt" \
,"Sample_167_5_K27Ac_Target_ChIPseeqer.txt" \
,"Sample_180_1_K27Ac_Target_ChIPseeqer.txt" \
,"Sample_180_2_K27Ac_Target_ChIPseeqer.txt" \
,"Sample_180_3_K27Ac_Target_ChIPseeqer.txt" \
,"Sample_180_4_K27Ac_Target_ChIPseeqer.txt" \
,"Sample_180_5_K27Ac_Target_ChIPseeqer.txt" \
,"Sample_181_4_K27Ac_Target_ChIPseeqer.txt"]

for i in file:
	for j in file:
		if i!=j:
			print "CompareIntervals -peakfile1 " +i+ " -peakfile2 "+j+" -show_ov_int 1 -outfile  "+i+j+'.overlap'
