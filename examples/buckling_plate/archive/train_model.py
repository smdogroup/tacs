import tensorflow as tf
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import niceplots

"""
GOAL: train a machine learning model on the buckling data as a surrogate model
Author: Sean Engelstad

notice the ML data in the ar-eig.png plot is similar to that of the "Theory of Elastic Stability"
lambda(AR) plot except that lambda does increase somewhat with higher AR. One major difference
is the Theory of Elastic Stability book uses Kirchoff-Love plate theory which excludes shear strain energy
and is still conservative. Whereas TACS includes shear strain energy in its geometrically exact plate theory soln.
"""

# get the dataset, first do a simple 1D dataset
df = pd.read_csv("buckle_data.csv")
AR = df["AR"].to_numpy()
eig1 = df["eig1"].to_numpy()
# print(AR)
# print(eig1)

# remove outliers, in this case any negative eigenvalues
x_train = []
y_train = []
for i in range(AR.shape[0]):
    if eig1[i] > 10.0:  # chosen min eigenvalue (dataset-dependent)
        x_train += [AR[i]]
        y_train += [eig1[i]]

# reshape into (N,1) arrays
x_train = np.array(x_train)
y_train = np.array(y_train)
n_train = x_train.shape[0]
x_train = np.reshape(x_train, newshape=(n_train, 1))
y_train = np.reshape(y_train, newshape=(n_train, 1))

# should split into training, validation, test data here
# but I'm not gonna do that just yet

# form the machine learning model
model = tf.keras.models.Sequential(
    [
        tf.keras.Input(shape=(1,)),
        tf.keras.layers.Dense(300, activation="tanh"),  #'relu'
        tf.keras.layers.Dense(200, activation="tanh"),  #'relu'
        tf.keras.layers.Dense(100, activation="tanh"),  #'relu'
        tf.keras.layers.Dense(1),
    ]
)
# tf.keras.layers.Dropout(0.2),

# https://github.com/keras-team/tf-keras/issues/501 - mean squared error worse in training step than in compiled version
optimizer = tf.keras.optimizers.Adam(learning_rate=1e-5)
loss_fn = tf.keras.losses.MeanSquaredError()
model.compile(
    optimizer=optimizer,
    loss=loss_fn,
)

history = model.fit(x_train, y_train, epochs=1000)

# add test + evaluation later
# model.evaluate(x_test,  y_test, verbose=2)

# plot the training history
# TODO : need to add validation set here
plt.style.use(niceplots.get_style())
plt.figure("training-hist")
plt.plot(history.history["loss"])
plt.title("training history")
plt.ylabel("loss")
plt.xlabel("epoch")
plt.yscale("log")
plt.savefig("training-hist.png", dpi=400)
plt.close("training-hist")

# plot the 1D model function
plt.figure("model-compare")
AR2 = np.linspace(1.0, 10.0, 300)
AR2 = np.reshape(AR2, newshape=(300, 1))
eig1_pred = model.predict(AR2)
plt.plot(x_train, y_train, label="data")
plt.plot(AR2, eig1_pred, label="model")
plt.legend()
plt.savefig("ar-eig-pred.png", dpi=400)
