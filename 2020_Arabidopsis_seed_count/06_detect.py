import sys,time
from pathlib import Path
import tensorflow as tf
flags = tf.app.flags
flags.DEFINE_string('base_path', '', 'the absolute path including graph_train directory')
flags.DEFINE_string('test_images', '', 'path of test images')
FLAGS = flags.FLAGS
from object_detection.utils import label_map_util
import numpy as np
import os
import six.moves.urllib as urllib
import tarfile
import tensorflow as tf
import zipfile
from collections import defaultdict
from io import StringIO
from PIL import Image
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ['CUDA_VISIBLE_DEVICES'] = "0"
from tensorflow.python.client import device_lib

config = tf.ConfigProto(log_device_placement=True)

#test images dir
PATH_TO_TEST_IMAGES_DIR = FLAGS.test_images
TEST_IMAGE_PATHS=[]
for root, dirs, files in os.walk(PATH_TO_TEST_IMAGES_DIR):
        for f in files:
            if not f.endswith("csv"):
                TEST_IMAGE_PATHS.append(os.path.join(PATH_TO_TEST_IMAGES_DIR,f))    
PATH_TO_CKPT=FLAGS.base_path+"/graph_train/frozen_inference_graph.pb"
PATH_TO_LABELS = FLAGS.base_path+"/mscoco_label_map.pbtxt"
NUM_CLASSES = 1
label_map = label_map_util.load_labelmap(PATH_TO_LABELS)
categories = label_map_util.convert_label_map_to_categories(label_map, max_num_classes=NUM_CLASSES, use_display_name=False)
category_index = label_map_util.create_category_index(categories)

def load_image_into_numpy_array(image):
  (im_width, im_height) = image.size
  return np.array(image.getdata()).reshape(
      (im_height, im_width, 3)).astype(np.uint8)

detection_graph = tf.Graph()
with detection_graph.as_default():
  od_graph_def = tf.GraphDef()
  with tf.gfile.GFile(PATH_TO_CKPT, 'rb') as fid:
    serialized_graph = fid.read()
    od_graph_def.ParseFromString(serialized_graph)
    tf.import_graph_def(od_graph_def, name='')

def run_inference_for_single_image(image, graph):
  with graph.as_default():
    with tf.Session(config = config) as sess:
      tensor_dict = {}
      tensor_name = 'detection_scores:0'
      tensor_dict["detection_scores"] = tf.get_default_graph().get_tensor_by_name(tensor_name)
      image_tensor = tf.get_default_graph().get_tensor_by_name('image_tensor:0')
      output_dict = sess.run(tensor_dict,feed_dict={image_tensor: np.expand_dims(image, 0)})
      output_dict['detection_scores'] = output_dict['detection_scores'][0]
  return output_dict
output = open(PATH_TO_TEST_IMAGES_DIR+"/Seed_count_results_"+time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())+".csv","a")
for image_path in TEST_IMAGE_PATHS:
  image = Image.open(image_path)
  image_np = load_image_into_numpy_array(image)
  image_np_expanded = np.expand_dims(image_np, axis=0)
  # Actual detection.
  output_dict = run_inference_for_single_image(image_np, detection_graph)
  scores = np.array(output_dict["detection_scores"])
  count = len(scores[scores>=0.5])
  output.write("%s,%d\n"%(image.filename.split("/")[-1],count))
output.close()
  
