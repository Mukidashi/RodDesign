work_dir=`pwd`

docker container run --runtime=nvidia  -it -w /work -v $work_dir/data:/work/data \
    -v $work_dir/src:/work/ \
    -v /tmp/.X11-unix/:/tmp/.X11-unix -e DISPLAY=$DISPLAY \
    --privileged --name rod_design rod_design_image bash
