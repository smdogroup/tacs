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
#print(AR)
#print(eig1)

# should split into training, validation, test data here
# but I'm not gonna do that just yet

# form the machine learning model
model = tf.keras.models.Sequential([
  tf.keras.Input(shape=(1,)),
  tf.keras.layers.Dense(128, activation='relu'),
  tf.keras.layers.Dropout(0.2),
  tf.keras.layers.Dense(20, activation='relu'),
  tf.keras.layers.Dense(1)
])
#tf.keras.layers.Dropout(0.2),

optimizer = tf.keras.optimizers.Adam(learning_rate=1e-4)
loss_fn = tf.keras.losses.MeanSquaredError()
model.compile(optimizer=optimizer,
              loss=loss_fn,)

model.fit(AR, eig1, epochs=300)

# add test + evaluation later
#model.evaluate(x_test,  y_test, verbose=2)

# plot the 1D model function
plt.style.use(niceplots.get_style())
AR2 = np.linspace(0.0, 10.0, 300)
eig1_pred = model.predict(AR2)
plt.plot(AR, eig1, label="data")
plt.plot(AR2, eig1_pred, label="model")
plt.savefig('ar-eig-pred.png', dpi=400)