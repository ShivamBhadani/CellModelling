# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def interpolate_waveform(csv_file, N, output_file=None):
    """
    Interpolates current vs. time data to create a waveform with N-times more time points.
    
    Parameters:
    - csv_file: Path to the input CSV file (should contain 'time' and 'current' columns).
    - N: Number of times more points to interpolate.
    - output_file: Path to save the interpolated data as a CSV (optional).
    
    Returns:
    - interpolated_df: DataFrame containing the interpolated time and current values.
    """
    # Load the data
    data = pd.read_csv(csv_file)
    time = data['time'].values
    current = data['Current'].values
    
    # Interpolation logic
    interp_func = interp1d(time, current, kind='linear')  # Linear interpolation
    new_time = np.linspace(time.min(), time.max(), N * len(time))
    new_current = interp_func(new_time)
    
    # Create a new DataFrame for interpolated data
    interpolated_df = pd.DataFrame({'time': new_time, 'Current': new_current})
    
    # Save interpolated data to CSV if specified
    if output_file:
        interpolated_df.to_csv(output_file, index=False)
    
    # Plot the original and interpolated waveform
    plt.figure(figsize=(10, 6))
    plt.plot(time, current, 'o', label='Original Data', markersize=4)
    plt.plot(new_time, new_current, '-', label='Interpolated Waveform')
    plt.xlabel('Time')
    plt.ylabel('Current')
    plt.title('Current vs. Time Waveform')
    plt.legend()
    plt.grid()
    plt.show()
    
    return interpolated_df

# Example usage
input_csv = 'la92shortdds.csv'  # Replace with your CSV file path
N = 10  # Number of times more points
output_csv = 'CurrentTimePerSecond.csv'  # Optional output file path

#interpolated_data = interpolate_waveform(input_csv, N, output_csv)

data = pd.read_csv(input_csv)
time = data['time'].values
current = data['Current'].values

# Interpolation logic
interp_func = interp1d(time, current, kind='linear')  # Linear interpolation
new_time = np.linspace(time.min(), time.max(), N * len(time))
new_current = interp_func(new_time)

# Create a new DataFrame for interpolated data
interpolated_df = pd.DataFrame({'time': new_time, 'Current': new_current})

# Save interpolated data to CSV if specified

interpolated_df.to_csv(output_csv, index=False)

# Plot the original and interpolated waveform
plt.figure(figsize=(10, 6))
plt.plot(time, current, 'o', label='Original Data', markersize=4)
plt.plot(new_time, new_current, '-', label='Interpolated Waveform')
plt.xlabel('Time')
plt.ylabel('Current')
plt.title('Current vs. Time Waveform')
plt.legend()
plt.grid()
plt.show()

print(interpolated_df)