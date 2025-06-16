#!/usr/bin/env python
"""
# Author: Jinzhao LI and Jishuai Miao
# File Name: __init__.py
# Description:
"""

__author__ = "Jinzhao LI and Jishuai Miao"
__email__ = "jinzhao9799@gmail.com"


from .Train_MultiGATE import train_MultiGATE
from .utils import Cal_Spatial_Net, Cal_gene_peak_Net, Cal_gene_peak_Net_new, Cal_gene_protein_Net, Stats_Spatial_Net, wnn_R, wnn_protein,mclust_R
from .model_MultiGATE import MGATE