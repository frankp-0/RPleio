s=$PWD"/../data/metain.txt.gz"
g=$PWD"/../data/sg.txt.gz"
e=$PWD"/../data/ce.txt.gz"
o=$PWD"/../results/pgc"
d=$PWD
mkdir ../results

# Run RPleio ====================
echo "Running RPleio"
Rscript ../R/pleio.R \
	-s $s \
	-g $g \
	-e $e \
	-n 1e5 \
	-c 10 \
	-o $o

echo "\n\n\n"

# Run PLEIO =====================
echo "Running PLEIO"
cd ~/pleio
./pleio \
    --metain $s \
    --sg $g \
    --ce $e \
    --nis 100000 \
    --parallel \
    --ncores 10 \
    --create \
    --out $o"_pleio"

echo "\n\n\n"

# Plot figures ==================
cd $d
echo "Plotting figures... \n"
Rscript "../R/results.R"
