import plotly.graph_objects as go

import numpy as np
import pandas as pd
import scipy

from scipy import signal

np.random.seed(1)

x = np.linspace(0, 10, 100)
y = np.sin(x)
noise = 2 * np.random.random(len(x)) - 1  # uniformly distributed between -1 and 1
y_noise = y + noise

y2 = signal.savgol_filter(y, 53,  # window size used for filtering
                         3)  # order of fitted polynomial

print(len(y), len(y2))
# fig = go.Figure()
# fig.add_trace(go.Scatter(
#     x=x,
#     y=y,
#     mode='markers',
#     marker=dict(size=2, color='black'),
#     name='Sine'
# ))

# fig.add_trace(go.Scatter(
#     x=x,
#     y=y_noise,
#     mode='markers',
#     marker=dict(
#         size=6,
#         color='royalblue',
#         symbol='circle-open'
#     ),
#     name='Noisy Sine'
# ))
#
# fig.add_trace(go.Scatter(
#     x=x,
#     y=signal.savgol_filter(y,
#                            53, # window size used for filtering
#                            3), # order of fitted polynomial
#     mode='markers',
#     marker=dict(
#         size=6,
#         color='mediumpurple',
#         symbol='triangle-up'
#     ),
#     name='Savitzky-Golay'
# ))
#
#
# # fig.show()
