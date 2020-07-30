# Arabidopsis seeds count using Tensorflow Faster-RCNN model

If you are interested in using Faster-RCNN to detect Arabidopsis seeds. Please see the following Github page:

https://github.com/FanruiMeng/Arabidopsis-Seed-Detection

# Training seed detection model

The following is for training new/updated Faster-RCNN model.

## 1. Tensorflow object detection API installation

* Tensorflow version 1.x (version 2 will not work)
  * [installation instruction](https://github.com/tensorflow/models/blob/master/research/object_detection/g3doc/installation.md)

## 2. Seed annotation

* This is only for training a new model. For applying the Faster R-CNN model, this is not necessary.
* We use [LabelImg](https://github.com/tzutalin/labelImg) to annotate our seeds. LabelImg generate a xml annotation file

<img src="https://github.com/FanruiMeng/Arabidopsis_seed_count/blob/master/Images/seeds_annotation.png?raw=true"  alt="Seed annotation" height="200"/>

* During the first round model training, we split one whole plate image into 4 quater images and annotate quater images mannually.
  
`python 00_split_scan_images.py image_directory_path`
  
## 3. Xml file transform to csv file

* Conversion script: 

`pyhton 02_xml_to_csv.py <need parameters>`
  
* The conversion results in a csv file:

  <table>
  <tr><td><b>filename</b></td> <td><b>width</b></td> <td><b>height</b></td> <td><b>class</b></td> <td><b>xmin</b></td><td><b>ymin</b></td><td><b>xmax</b></td><td><b>ymax</b></td></tr>
  <tr><td>scan21_111918022.jpg</td> <td>2900</td> <td>2900</td> <td>seed</td> <td>1325</td><td>813</td><td>1352</td><td>837</td></tr>
  <tr><td>scan21_111918022.jpg</td> <td>2900</td> <td>2900</td> <td>seed</td> <td>664</td><td>1094</td><td>691</td><td>1116</td></tr>
  <tr><td>..</td> <td>..</td> <td>..</td> <td>..</td> <td>..</td><td>..</td><td>..</td><td>..</td></tr>
  </table>
  
## 4. Convert CSV file to Tensorflow tfrecord file

* Conversion script: 

`python 02_generate_tfrecord.py --csv_input=annotation/seeds_labels.csv --output_path=train.record`

## 5. Download tensorflow object detection API pre-trained Faster R-CNN model

* In your terminal, issue the following command:

`wget http://download.tensorflow.org/models/object_detection/faster_rcnn_inception_v2_coco_2018_01_28.tar.gz`

`unzip faster_rcnn_inception_v2_coco_2018_01_28.tar.gz`

## 6. Pipeline configuration

### Input configuration

* In the unzip folder `faster_rcnn_inception_v2_coco_2018_01_28`
* Replace the `pipeline.config` file with [the one provided in this repository](https://github.com/ShiuLab/Manuscript_Code/blob/master/2020_Arabidopsis_seed_count/pipeline.config)
* For more information on how to change the configuration, see [this document](https://github.com/tensorflow/models/blob/master/research/object_detection/g3doc/configuring_jobs.md).

## 7. Model training

* Train model with the following parameters
  * `logtostderr`: provide logs during training
  * `pipeline_config_path`: the path to the configuration file
  * `train_dir`: output directory for the model
  * `num_clones`: number of processing units used to train the model

`python3 03_train.py --logtostderr --pipeline_config_path=pipeline.config --train_dir=train_dir --num_clones=3`

## 8. Generate a frozen (final) model

* Geneate with the following parameters:
  * `input_type`: the type of input
  * `pipeline_config_path`: the path to the configuration file
  * `trained_checkpoint_prefix`: prefix of the names for the model checkpoint files
  * `output_directory`: name for output directory

`python3 05_export_inference_graph.py --input_type image_tensor --pipeline_config_path pipeline.config --trained_checkpoint_prefix train_dir/model.ckpt- --output_directory graph_train`

## 9. Detect seeds using trained model

* Parameter: <need info>

`python 06_detect.py <need info>`

## 10. Accuracy measurement

* Measure accuracy, precision, recall and f1 at IOU 0.5 using 07_01_accuracy_measurement.py <br>

`python 07_01_accuracy_measurement.py ground.csv detected.csv`

## 11. Estimating seed density

* Determine the average seed number in a circle with a radius of 30 pixels.

`Rscript seed_density.r`
