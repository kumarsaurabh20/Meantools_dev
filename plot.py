import pandas as pd
import matplotlib.pyplot as plt

def frequency_distribution(counts, kind, color, rotate=0, ha='', xlabel='', ylabel='', title=''):
	# Plotting the frequency distribution as a bar plot
	plt.figure(figsize=(10,6))
	value_counts.plot(kind=kind, color=color)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)
	# Rotate x-axis labels for better readability
	plt.xticks(rotation=rotate, ha=ha)
	plt.tight_layout()  # Adjust layout to prevent clipping of labels
	plt.show()


if __name__ == "__main__":

	df = pd.read_csv("dummy_transition.csv")
	value_counts = df['mass_transition'].value_counts().sort_index()
	frequency_distribution(value_counts, 'bar', 'blue', rotate=45, ha='right', xlabel='Mass Transition Value', ylabel='Frequency', title='Distribution of Mass Transitions')
