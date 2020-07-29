<h2>Arabidopsis seeds count using Tensorflow Faster-RCNN model</h2>
<h3>1 Tensorflow object detection API installation</h3>
  tensorflow version lower than 2.0.
  please see: 
  https://github.com/tensorflow/models/blob/master/research/object_detection/g3doc/installation.md

<h3>2 Seed annotation</h3>
  we use LabelImg annotate our seeds, LabelImg generate a xml annotation file.<br>
  LabelImg please see: https://github.com/tzutalin/labelImg<br>
  For our project, we split one whole plate image into 4 quater images and annotate quater images mannually.<br>
  <b>a. split image:</b><br>
<i>python 00_split_scan_images.py</i><br>
<b>b. seed annotation</b>
  <img src="https://github.com/FanruiMeng/Arabidopsis_seed_count/blob/master/Images/seeds_annotation.png?raw=true"  alt="Seed annotation" height="200" width="300"/>
<h3>3 Xml file transform to csv file</h3>
  <i>pyhton 02_xml_to_csv.py</i><br>
  
  Results is a csv file, like this:<br><br>
  <table>
  <tr><td><b>filename</b></td> <td><b>width</b></td> <td><b>height</b></td> <td><b>class</b></td> <td><b>xmin</b></td><td><b>ymin</b></td><td><b>xmax</b></td><td><b>ymax</b></td></tr>
  <tr><td>scan21_111918022.jpg</td> <td>2900</td> <td>2900</td> <td>seed</td> <td>1325</td><td>813</td><td>1352</td><td>837</td></tr>
  <tr><td>scan21_111918022.jpg</td> <td>2900</td> <td>2900</td> <td>seed</td> <td>664</td><td>1094</td><td>691</td><td>1116</td></tr>
  <tr><td>..</td> <td>..</td> <td>..</td> <td>..</td> <td>..</td><td>..</td><td>..</td><td>..</td></tr>
  </table>
<h3>4. seed_labels.csv transform to tensorflow tfrecord file </h3>

  <i>python 02_generate_tfrecord.py --csv_input=annotation/seeds_labels.csv --output_path=train.record</i>

<h3>5 Download tensorflow object detection api pre-trained faster rcnn model to work directory.</h3>
  wget http://download.tensorflow.org/models/object_detection/faster_rcnn_inception_v2_coco_2018_01_28.tar.gz<br>
  unzip faster_rcnn_inception_v2_coco_2018_01_28.tar.gz

<h3>6 pipeline configuration
  <h4>a. input configuration</h4>
  train_input_reader: {<br>
  &nbsp;&nbsp;tf_record_input_reader {<br>
   &nbsp;&nbsp;&nbsp;&nbsp;input_path: "train.record"<br>
    &nbsp;&nbsp;}<br>
  
  &nbsp;&nbsp;label_map_path: "mscoco_label_map.pb<br>
}
  <h4>b. The label_map.pbtxt file like below:</h4>
  item {<br>
    &nbsp;&nbsp;id: 1<br>
    &nbsp;&nbsp;name: "seed"<br>
  }<br>
  c. Another configurations please see: 
  https://github.com/tensorflow/models/blob/master/research/object_detection/g3doc/configuring_jobs.md
<h3>7 Model training</h3>
  <i>python3 03_train.py --logtostderr --pipeline_config_path=pipeline.config --train_dir=train_dir --num_clones=3</i>
<h3>8 Generate frozen model </h3>
  <i>python3 05_export_inference_graph.py --input_type image_tensor --pipeline_config_path pipeline.config --trained_checkpoint_prefix train_dir/model.ckpt- --output_directory graph_train</i>
<h3>9 Detect seeds using trained model </h3>
  <i>python 06_detect.py</i>
<h3>10 Accuracy measurement</h3>
Measure accuracy, precision, recall and f1 at IOU 0.5 using 07_01_accuracy_measurement.py <br>
<i>python 07_01_accuracy_measurement.py ground.csv detected.csv</i>
<h3>11 seed density</h3>
Average seed number in a circle with a radius of 30 pixels.<br>
<i>Rscript seed_density.r</i>
