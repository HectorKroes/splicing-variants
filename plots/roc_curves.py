import numpy as np
import matplotlib.pyplot as plt

# Sensitivity and specificity values for each curve
SenSpliceAI = [1, 0.86, 0.842, 0.829, 0.817, 0.809, 0.801, 0.794, 0.785, 0.776, 0.768, 0.756, 0.745, 0.732, 0.716, 0.697, 0.673, 0.642, 0.596, 0.509, 0]
SpeSpliceAI = [0, 0.888, 0.942, 0.964, 0.974, 0.981, 0.985, 0.988, 0.99, 0.992, 0.992, 0.993, 0.994, 0.995, 0.995, 0.996, 0.996, 0.997, 0.997, 0.998, 1]

SenSquirls = [1, 0.767, 0.751, 0.744, 0.74, 0.737, 0.732, 0.728, 0.723, 0.718, 0.711, 0.704, 0.695, 0.683, 0.677, 0.67, 0.664, 0.657, 0.645, 0.564, 0]
SpeSquirls = [0, 0.953, 0.969, 0.973, 0.976, 0.979, 0.981, 0.982, 0.984, 0.985, 0.986, 0.987, 0.988, 0.989, 0.99, 0.99, 0.991, 0.992, 0.993, 0.994, 1]

SenPipeline = [1.0, 0.864, 0.857, 0.856, 0.854, 0.854, 0.837, 0.836, 0.835, 0.834, 0.833, 0.818, 0.816, 0.814, 0.812, 0.81, 0.784, 0.781, 0.767, 0.73, 0]
SpePipeline = [0.0, 0.936, 0.95, 0.954, 0.957, 0.959, 0.975, 0.977, 0.978, 0.979, 0.98, 0.984, 0.985, 0.986, 0.987, 0.988, 0.99, 0.991, 0.991, 0.993, 1]

# Calculate false positive rate and true positive rate from sensitivity and specificity
fprSpliceAI = [1 - spec for spec in SpeSpliceAI]
tprSpliceAI = SenSpliceAI

fprSquirls = [1 - spec for spec in SpeSquirls]
tprSquirls = SenSquirls

fprPipeline = [1 - spec for spec in SpePipeline]
tprPipeline = SenPipeline

# Sort the points in ascending order of false positive rate
pointsSpliceAI = sorted(zip(fprSpliceAI, tprSpliceAI))
pointsSquirls = sorted(zip(fprSquirls, tprSquirls))
pointsPipeline = sorted(zip(fprPipeline, tprPipeline))

# Calculate the ROC curve and AUC using the trapezoidal rule
rocSpliceAI = np.trapz([tpr for _, tpr in pointsSpliceAI], [fpr for fpr, _ in pointsSpliceAI])
rocSquirls = np.trapz([tpr for _, tpr in pointsSquirls], [fpr for fpr, _ in pointsSquirls])
rocPipeline = np.trapz([tpr for _, tpr in pointsPipeline], [fpr for fpr, _ in pointsPipeline])

# Plot the ROC curves
plt.figure()
plt.plot(fprPipeline, tprPipeline, label='Our Pipeline (AUC = %0.3f)' % rocPipeline)
plt.plot(fprSpliceAI, tprSpliceAI, label='SpliceAI (AUC = %0.3f)' % rocSpliceAI)
plt.plot(fprSquirls, tprSquirls, label='SQUIRLS (AUC = %0.3f)' % rocSquirls)
plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Benchmarking of splice-altering variant prediction tools')
plt.legend(loc="lower right")
plt.show()