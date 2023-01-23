
workdir=$1
webdir=$2

echo "wordir=" $workdir
echo "webdir=" $webdir

cd $workdir
# tar files 
tar cvf view.tar viewer figure

# copy to webdir
cp view.tar ${webdir}

# untar at webdir
cd ${webdir}
tar xvf view.tar

# change ownership to 755
chmod -R 755 viewer
chmod -R 755 figure 

# remove view.tar
rm -f view.tar 
