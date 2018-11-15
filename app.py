from flask import Flask, render_template, request, redirect
#from six.moves import urllib
import numpy as np
import matplotlib.pyplot as plt
import os

import requests
from collections import defaultdict


from bokeh.embed import components
from bokeh.plotting import figure
from bokeh.resources import INLINE
from bokeh.util.string import encode_utf8
from bokeh.models.graphs import from_networkx

import networkx as nx

from bokeh.layouts import column
from bokeh.models import *
from bokeh.plotting import Figure, output_file, show
from bokeh.io import output_notebook, reset_output
from bokeh.io import show, output_file
from bokeh.plotting import figure
from bokeh.palettes import Spectral4


import gzip, glob, os, scipy
import pickle as pick
import numpy as np
def pickle(fname, obj):
    pick.dump(obj=obj, file=gzip.open(fname, "wb", compresslevel=3), protocol=2)
def unpickle(fname):
    return pick.load(gzip.open(fname, "rb"))


disease_network=unpickle('disease_network.obj')
drug_affinity_dict=unpickle('drug_affinity_dict.obj')

all_disease = ['eye disease', 'bladder disease', 'skeletal system disease', 'reproductive system disease', 'immune system disease', 'endocrine system disease', 'genetic disorder',
 'hematological system disease', 'skin disease', 'biological process', 'neoplasm', 'cardiovascular disease', 'nervous system disease', 'metabolic disease', 'respiratory system disease', 'head and neck disorder', 'infectious disease', 'kidney disease', 'other disease', 'digestive system disease']

uni_list = ['P03372','P24941','P00734']


def construct_prot_disease_network(disease_list,protein_list,disease_network):
  orig_nodes = set([s for s,_ in disease_network.nodes.items()])
  disease_list = list(set(disease_list))
  protein_list = list(set(protein_list))
  sub_net = nx.Graph()
  add_list = set()
  neighbor_list = set()
  for p in list(orig_nodes.intersection(protein_list)):
      add_list.add(p)
      neighbor_list.add(p)
  for d in list(orig_nodes.intersection(disease_list)):
      add_list.add(d)
      neighbor_list.add(d)
  for node in add_list:
      for nn in list(disease_network.neighbors(node)):
          neighbor_list.add(nn)
  for node in list(neighbor_list):
      sub_net.add_node(node)
  
  for node in list(neighbor_list):
      for nnode in list(disease_network.neighbors(node)):
          if not sub_net.has_edge(node,nnode):
              sub_net.add_edge(node,nnode,weight=disease_network.get_edge_data(node,nnode)['weight'])
  elarge = [(u, v) for (u, v, d) in sub_net.edges(data=True) if d['weight'] > 0.1]
  esmall = [(u, v) for (u, v, d) in sub_net.edges(data=True) if d['weight'] <= 0.1]
  pos_map = nx.spring_layout(sub_net)
  color_map=[[],[]]
  prot_list,dis_list=[],[]
  label_prot={}
  for node in sub_net.nodes():
      if node in uni_list:
          color_map[0].append('blue')
          prot_list.append(node)
          label_prot[node]=node
      elif node in all_disease:
          color_map[1].append('green')
          disease_list.append(node)

  return [sub_net,pos_map,color_map,prot_list,disease_list,esmall,elarge,label_prot]

def construct_prot_drug_network(protein_list,drug_affinity_dict,top_n):
  all_drug_list=set()
  min_drug_set={}
  for p in protein_list:
      drug_list=[]
      min_drug_set[p]=[]
      for d in drug_affinity_dict[p].keys():
          drug_list.append((drug_affinity_dict[p][d],d))
      drug_list=sorted(drug_list,reverse=True)
      if len(drug_list)>top_n:
          for d in drug_list[0:top_n]:
              all_drug_list.add(d[1])
              min_drug_set[p].append(d)
      else:
          for d in drug_list:
              all_drug_list.add(d)      
              min_drug_set[p].append(d)
        
  sub_net = nx.Graph()
  for p in protein_list:
      sub_net.add_node(p)
  for d in list(all_drug_list):
      sub_net.add_node(d)
  for p in protein_list:
      for aff,d in min_drug_set[p]:
          sub_net.add_edge(p,d,weight=aff)
  return sub_net

app = Flask(__name__)

app.vars={}


@app.route('/')
def index():
  return render_template('index.html')

@app.route('/about')
def about():
  return render_template('about.html')

@app.route('/disease_select',methods=['GET','POST'])
def sel_disease():
  return render_template('disease_select.html')

@app.route('/drug_select',methods=['GET','POST'])
def sel_drug():
  return render_template('drug_select.html')

@app.route('/graph_data_prot',methods=['GET','POST'])
def make_graph_prot():
    
  protein_list = [p for p in request.form['protein_list'].split(',') if len(p)>0]
  disease_list = [p for p in request.form['disease_list'].split(',') if len(p)>0]
    
  sub_net,pos_map,color_map,prot_list,disease_list,esmall,elarge,label_prot=construct_prot_disease_network(disease_list,protein_list,disease_network)
  
    
  fig = figure(plot_width=400, plot_height=400,x_range=Range1d(-1.1,1.1), y_range=Range1d(-1.1,1.1))

    
  graph_renderer = from_networkx(sub_net, nx.spring_layout, scale=1, center=(0,0))
    
  edge_color_list=[]
  for e1,e2 in zip(graph_renderer.edge_renderer.data_source.data['start'],graph_renderer.edge_renderer.data_source.data['end']):
    if sub_net.get_edge_data(e1,e2)['weight']>0.1:
      edge_color_list.append('#000000')
    else:
      edge_color_list.append('#ff4c4c')
  color_list=[]
  size_list=[]
  for i in graph_renderer.node_renderer.data_source.data['index']:
    if i in uni_list:
      color_list.append('#3232ff')
      size_list.append(20)
    elif i in all_disease:
      color_list.append('#228B22')
      size_list.append(15)
  graph_renderer.node_renderer.data_source.data['color']=color_list

#  time_b =  request.form['time_b']    

  fig.title.text = "Graph Interaction Demonstration"

  hover = HoverTool(tooltips=[("", "@name")])
  graph_renderer.node_renderer.data_source.data['name']=graph_renderer.node_renderer.data_source.data['index']  
  graph_renderer.node_renderer.data_source.data['ncolor']=color_list
  graph_renderer.node_renderer.data_source.data['nsize']=size_list
 
  graph_renderer.node_renderer.glyph = Circle(size='nsize', fill_color='ncolor')
  graph_renderer.node_renderer.selection_glyph = Circle(size='nsize', fill_color=Spectral4[2])
  graph_renderer.node_renderer.hover_glyph = Circle(size='nsize', fill_color='ncolor')
  graph_renderer.edge_renderer.data_source.data['lcolor']=edge_color_list
  graph_renderer.edge_renderer.glyph = MultiLine(line_color='lcolor', line_alpha=0.6, line_width=3)
  graph_renderer.edge_renderer.selection_glyph = MultiLine(line_color=Spectral4[2], line_width=3)
  graph_renderer.edge_renderer.hover_glyph = MultiLine(line_color='lcolor', line_width=3)

  graph_renderer.selection_policy = NodesAndLinkedEdges()
  graph_renderer.inspection_policy = EdgesAndLinkedNodes()
  fig.add_tools(hover, TapTool(), BoxSelectTool(), WheelZoomTool(), LassoSelectTool(),PanTool())
  graph_renderer.inspection_policy = NodesAndLinkedEdges()

  fig.renderers.append(graph_renderer)
  js_resources = INLINE.render_js()
  css_resources = INLINE.render_css()
  script, div = components(fig)
  html = render_template('graph_data_prot.html',plot_script=script,plot_div=div,js_resources=js_resources,css_resources=css_resources,)
  return encode_utf8(html)
    
    
    
@app.route('/graph_data_drug',methods=['GET','POST'])
def make_graph_drug():
    
  protein_list = [p for p in request.form['protein_list'].split(',') if len(p)>0]
  top_n = int(request.form['num_drugs'])

  sub_net=construct_prot_drug_network(protein_list,drug_affinity_dict,top_n)
  
  fig = figure(plot_width=400, plot_height=400,x_range=Range1d(-1.1,1.1), y_range=Range1d(-1.1,1.1))

  graph_renderer = from_networkx(sub_net, nx.spring_layout, scale=1, center=(0,0))
    
  edge_color_list=[]
  for e1,e2 in zip(graph_renderer.edge_renderer.data_source.data['start'],graph_renderer.edge_renderer.data_source.data['end']):
    if sub_net.get_edge_data(e1,e2)['weight']>8.0:
      edge_color_list.append('#000000')
    else:
      edge_color_list.append('#ff4c4c')
  color_list=[]
  size_list=[]
  for i in graph_renderer.node_renderer.data_source.data['index']:
    if i in uni_list:
      color_list.append('#3232ff')
      size_list.append(20)
    else:
      color_list.append('#228B22')
      size_list.append(10)
  graph_renderer.node_renderer.data_source.data['color']=color_list

  fig.title.text = "Graph Interaction Demonstration"

  hover = HoverTool(tooltips=[("", "@name")])
  graph_renderer.node_renderer.data_source.data['name']=graph_renderer.node_renderer.data_source.data['index']  
  graph_renderer.node_renderer.data_source.data['ncolor']=color_list
  graph_renderer.node_renderer.data_source.data['nsize']=size_list
 
  graph_renderer.node_renderer.glyph = Circle(size='nsize', fill_color='ncolor')
  graph_renderer.node_renderer.selection_glyph = Circle(size='nsize', fill_color=Spectral4[2])
  graph_renderer.node_renderer.hover_glyph = Circle(size='nsize', fill_color='ncolor')
  graph_renderer.edge_renderer.data_source.data['lcolor']=edge_color_list
  graph_renderer.edge_renderer.glyph = MultiLine(line_color='lcolor', line_alpha=0.6, line_width=3)
  graph_renderer.edge_renderer.selection_glyph = MultiLine(line_color=Spectral4[2], line_width=3)
  graph_renderer.edge_renderer.hover_glyph = MultiLine(line_color='lcolor', line_width=3)

  graph_renderer.selection_policy = NodesAndLinkedEdges()
  graph_renderer.inspection_policy = EdgesAndLinkedNodes()
  fig.add_tools(hover, TapTool(), BoxSelectTool(), WheelZoomTool(), LassoSelectTool(),PanTool())
  graph_renderer.inspection_policy = NodesAndLinkedEdges()

  fig.renderers.append(graph_renderer)
  js_resources = INLINE.render_js()
  css_resources = INLINE.render_css()
  script, div = components(fig)
  html = render_template('graph_data_drug.html',plot_script=script,plot_div=div,js_resources=js_resources,css_resources=css_resources,)
  return encode_utf8(html)
    

if __name__ == '__main__':
#  app.run(port=5000,debug=True)
  app.run(port=33507)
