grep -w -F -f OFF_MEF_IVF_NT_1C.csv TSS.bed > OFF_MEF_IVF_NT_1C.bed
grep -w -F -f OFF_MEF_IVF_NT_2C.csv TSS.bed > OFF_MEF_IVF_NT_2C.bed
grep -w -F -f OFF_MEF_IVF_NT_4C.csv TSS.bed > OFF_MEF_IVF_NT_4C.bed
grep -w -F -f OFF_MEF_IVF_NT_8C.csv TSS.bed > OFF_MEF_IVF_NT_8C.bed
grep -w -F -f OFF_MEF_IVF_NT_B.csv TSS.bed > OFF_MEF_IVF_NT_B.bed
grep -w -F -f OFF_MEF_IVF_NT_M.csv TSS.bed > OFF_MEF_IVF_NT_M.bed
grep -w -F -f ON_MEF_IVF_NT_1C.csv TSS.bed > ON_MEF_IVF_NT_1C.bed
grep -w -F -f ON_MEF_IVF_NT_2C.csv TSS.bed > ON_MEF_IVF_NT_2C.bed
grep -w -F -f ON_MEF_IVF_NT_4C.csv TSS.bed > ON_MEF_IVF_NT_4C.bed
grep -w -F -f ON_MEF_IVF_NT_8C.csv TSS.bed > ON_MEF_IVF_NT_8C.bed
grep -w -F -f ON_MEF_IVF_NT_B.csv TSS.bed > ON_MEF_IVF_NT_B.bed
grep -w -F -f ON_MEF_IVF_NT_M.csv TSS.bed > ON_MEF_IVF_NT_M.bed
grep -w -F -f RD_MEF_IVF_NT_1C.csv TSS.bed > RD_MEF_IVF_NT_1C.bed
grep -w -F -f RD_MEF_IVF_NT_2C.csv TSS.bed > RD_MEF_IVF_NT_2C.bed
grep -w -F -f RD_MEF_IVF_NT_4C.csv TSS.bed > RD_MEF_IVF_NT_4C.bed
grep -w -F -f RD_MEF_IVF_NT_8C.csv TSS.bed > RD_MEF_IVF_NT_8C.bed
grep -w -F -f RD_MEF_IVF_NT_B.csv TSS.bed > RD_MEF_IVF_NT_B.bed
grep -w -F -f RD_MEF_IVF_NT_M.csv TSS.bed > RD_MEF_IVF_NT_M.bed
grep -w -F -f RU_MEF_IVF_NT_1C.csv TSS.bed > RU_MEF_IVF_NT_1C.bed
grep -w -F -f RU_MEF_IVF_NT_2C.csv TSS.bed >  RU_MEF_IVF_NT_2C.bed
grep -w -F -f RU_MEF_IVF_NT_4C.csv TSS.bed > RU_MEF_IVF_NT_4C.bed
grep -w -F -f RU_MEF_IVF_NT_8C.csv TSS.bed > RU_MEF_IVF_NT_8C.bed
grep -w -F -f RU_MEF_IVF_NT_B.csv TSS.bed > RU_MEF_IVF_NT_B.bed
grep -w -F -f RU_MEF_IVF_NT_M.csv TSS.bed > RU_MEF_IVF_NT_M.bed

awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' OFF_MEF_IVF_NT_1C.bed > OFF_MEF_IVF_NT_1C.1.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' OFF_MEF_IVF_NT_2C.bed > OFF_MEF_IVF_NT_2C.1.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' OFF_MEF_IVF_NT_4C.bed > OFF_MEF_IVF_NT_4C.1.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' OFF_MEF_IVF_NT_8C.bed > OFF_MEF_IVF_NT_8C.1.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' OFF_MEF_IVF_NT_B.bed > OFF_MEF_IVF_NT_B.1.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' OFF_MEF_IVF_NT_M.bed > OFF_MEF_IVF_NT_M.1.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' ON_MEF_IVF_NT_1C.bed > ON_MEF_IVF_NT_1C.t.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' ON_MEF_IVF_NT_2C.bed > ON_MEF_IVF_NT_2C.1.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' ON_MEF_IVF_NT_4C.bed >ON_MEF_IVF_NT_4C.1.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' ON_MEF_IVF_NT_8C.bed >ON_MEF_IVF_NT_8C.1.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' ON_MEF_IVF_NT_B.bed >ON_MEF_IVF_NT_B.1.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' ON_MEF_IVF_NT_M.bed > ON_MEF_IVF_NT_M.1.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' RD_MEF_IVF_NT_1C.bed >RD_MEF_IVF_NT_1C.1.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' RD_MEF_IVF_NT_2C.bed >RD_MEF_IVF_NT_2C.1.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' RD_MEF_IVF_NT_4C.bed > RD_MEF_IVF_NT_4C.1.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' RD_MEF_IVF_NT_8C.bed > RD_MEF_IVF_NT_8C.1.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' RD_MEF_IVF_NT_B.bed > RD_MEF_IVF_NT_B.1.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' RD_MEF_IVF_NT_M.bed > RD_MEF_IVF_NT_M.1.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' RU_MEF_IVF_NT_1C.bed > RU_MEF_IVF_NT_1C.1.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' RU_MEF_IVF_NT_2C.bed > RU_MEF_IVF_NT_2C.1.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' RU_MEF_IVF_NT_4C.bed > RU_MEF_IVF_NT_4C.1.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' RU_MEF_IVF_NT_8C.bed > RU_MEF_IVF_NT_8C.1.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' RU_MEF_IVF_NT_B.bed > RU_MEF_IVF_NT_B.1.bed
awk '{ print $2 "\t" $3 "\t" $4 "\t" $1}' RU_MEF_IVF_NT_M.bed > RU_MEF_IVF_NT_M.1.bed


bedtools slop -i OFF_MEF_IVF_NT_1C.1.bed -g chr.sizes -b 1000 > OFF_MEF_IVF_NT_1C.t.bed
bedtools slop -i OFF_MEF_IVF_NT_2C.1.bed -g chr.sizes -b 1000 >  OFF_MEF_IVF_NT_2C.t.bed
bedtools slop -i OFF_MEF_IVF_NT_4C.1.bed -g chr.sizes -b 1000 > OFF_MEF_IVF_NT_4C.t.bed
bedtools slop -i OFF_MEF_IVF_NT_8C.1.bed -g chr.sizes -b 1000 > OFF_MEF_IVF_NT_8C.t.bed
bedtools slop -i OFF_MEF_IVF_NT_B.1.bed -g chr.sizes -b 1000 > OFF_MEF_IVF_NT_B.t.bed
bedtools slop -i OFF_MEF_IVF_NT_M.1.bed -g chr.sizes -b 1000 > OFF_MEF_IVF_NT_M.t.bed
bedtools slop -i ON_MEF_IVF_NT_1C.1.bed -g chr.sizes -b 1000 > ON_MEF_IVF_NT_1C.t.bec
bedtools slop -i ON_MEF_IVF_NT_2C.1.bed -g chr.sizes -b 1000 > ON_MEF_IVF_NT_2C.t.bed
bedtools slop -i ON_MEF_IVF_NT_4C.1.bed -g chr.sizes -b 1000 > ON_MEF_IVF_NT_4C.t.bed
bedtools slop -i ON_MEF_IVF_NT_8C.1.bed -g chr.sizes -b 1000 > ON_MEF_IVF_NT_8C.t.bed
bedtools slop -i ON_MEF_IVF_NT_B.1.bed -g chr.sizes -b 1000 > ON_MEF_IVF_NT_B.t.bed
bedtools slop -i ON_MEF_IVF_NT_M.1.bed -g chr.sizes -b 1000 > ON_MEF_IVF_NT_M.t.bed
bedtools slop -i RD_MEF_IVF_NT_1C.1.bed -g chr.sizes -b 1000 > RD_MEF_IVF_NT_1C.t.bed
bedtools slop -i RD_MEF_IVF_NT_2C.1.bed -g chr.sizes -b 1000 > RD_MEF_IVF_NT_2C.t.bed
bedtools slop -i RD_MEF_IVF_NT_4C.1.bed -g chr.sizes -b 1000 > RD_MEF_IVF_NT_4C.t.bed
bedtools slop -i RD_MEF_IVF_NT_8C.1.bed -g chr.sizes -b 1000 > RD_MEF_IVF_NT_8C.t.bed
bedtools slop -i RD_MEF_IVF_NT_B.1.bed -g chr.sizes -b 1000 > RD_MEF_IVF_NT_B.t.bed
bedtools slop -i RD_MEF_IVF_NT_M.1.bed -g chr.sizes -b 1000 > RD_MEF_IVF_NT_M.t.bed
bedtools slop -i RU_MEF_IVF_NT_1C.1.bed -g chr.sizes -b 1000 > RU_MEF_IVF_NT_1C.t.bed
bedtools slop -i RU_MEF_IVF_NT_2C.1.bed -g chr.sizes -b 1000 > RU_MEF_IVF_NT_2C.t.bed
bedtools slop -i RU_MEF_IVF_NT_4C.1.bed -g chr.sizes -b 1000 > RU_MEF_IVF_NT_4C.t.bed
bedtools slop -i RU_MEF_IVF_NT_8C.1.bed -g chr.sizes -b 1000 > RU_MEF_IVF_NT_8C.t.bed
bedtools slop -i RU_MEF_IVF_NT_B.1.bed -g chr.sizes -b 1000 > RU_MEF_IVF_NT_B.t.bed
bedtools slop -i RU_MEF_IVF_NT_M.1.bed -g chr.sizes -b 1000 > RU_MEF_IVF_NT_M.t.bed

