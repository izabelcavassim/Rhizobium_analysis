qx --no-scratch -t 00:30:00 -c 16 -w "/home/agb/tools/proteinortho_v5.13/proteinortho5.pl -cpus=16 -synteny -step=1 data/*.faa"

for i in {1..199}; do
    qx --no-scratch -t 24:00:00 -c 4 "/home/agb/tools/proteinortho_v5.13/proteinortho5.pl -cpus=4 -synteny -step=2 -startat=$i -stopat=$i data/*.faa"
done

qx --no-scratch -t 4:00:00 -c 8 "/home/agb/tools/proteinortho_v5.13/proteinortho5.pl -cpus=8 -synteny -step=3 data/*.faa"