
import streamlit as st
import pandas as pd
import os
from venny4py.venny4py import venny4py

# Sample placeholder sets and keys
keys = ["File1", "File2", "File3", "File4"]
sets = [
    {"A", "B", "C"},
    {"B", "C", "D"},
    {"C", "D", "E"},
    {"A", "C", "E"}
]

# Corrected function call with required 'labels' parameter
plot = venny4py(labels=keys, sets=sets)
st.pyplot(plot.gcf())
