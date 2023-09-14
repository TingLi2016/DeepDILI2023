param($image, $sdfFile)

docker run -i --rm $image /deep-dili/main_cli.py /data/training/online-objs.pkl $sdfFile
