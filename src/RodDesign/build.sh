cd /work
rm -rf RodDesign_build
mkdir RodDesign_build

cd RodDesign_build
cmake -DCMAKE_BUILD_TYPE=Debug ../RodDesign && make 
mv ../RodDesign_build/RodDesignGUI ../RodDesign