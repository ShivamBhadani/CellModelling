# -*- coding: utf-8 -*-
"""
Created on Sat Mar 29 20:47:38 2025

@author: iitm9
"""

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

class LookupTableEstimator:
    def __init__(self, file_path):
        # Load the lookup table once during initialization
        self.lookup_table = pd.read_csv(file_path)
        
        # Extract column headers
        headers = self.lookup_table.columns.tolist()
        
        # First column: input values
        self.inputs = self.lookup_table[headers[0]].values
        
        # Remaining columns: output data
        self.temperatures = [int(header.split(" ")[-1][:-1]) for header in headers[1:]]
        self.data = self.lookup_table.iloc[:, 1:].values

    def output(self, input_value, temperature):
        # Interpolate/extrapolate for the input value across all temperatures
        if input_value < self.inputs.min() or input_value > self.inputs.max():
            # Extrapolation for input if outside bounds
            input_func = interp1d(self.inputs, self.data, axis=0, fill_value="extrapolate")
            outputs_at_temperatures = input_func(input_value)
        else:
            # Interpolation for input if within bounds
            input_func = interp1d(self.inputs, self.data, axis=0)
            outputs_at_temperatures = input_func(input_value)
        
        # Interpolate/extrapolate for the temperature
        if temperature < self.temperatures[0] or temperature > self.temperatures[-1]:
            # Extrapolation for temperature if outside bounds
            temp_func = interp1d(self.temperatures, outputs_at_temperatures, fill_value="extrapolate")
            output = temp_func(temperature)
        else:
            # Interpolation for temperature if within bounds
            temp_func = interp1d(self.temperatures, outputs_at_temperatures)
            output = temp_func(temperature)
        
        return output

    def estimate_input(self, output_value, temperature):
        # Interpolate/extrapolate across temperatures to get outputs at all inputs
        if temperature < self.temperatures[0] or temperature > self.temperatures[-1]:
            # Extrapolation for temperature if outside bounds
            temp_func = interp1d(self.temperatures, self.data.T, axis=0, fill_value="extrapolate")
            outputs_at_inputs = temp_func(temperature)
        else:
            # Interpolation for temperature if within bounds
            temp_func = interp1d(self.temperatures, self.data.T, axis=0)
            outputs_at_inputs = temp_func(temperature)
        
        # Interpolate/extrapolate for the input value
        if output_value < outputs_at_inputs.min() or output_value > outputs_at_inputs.max():
            # Extrapolation for input if outside bounds
            input_func = interp1d(outputs_at_inputs, self.inputs, fill_value="extrapolate")
            input_value = input_func(output_value)
        else:
            # Interpolation for input if within bounds
            input_func = interp1d(outputs_at_inputs, self.inputs)
            input_value = input_func(output_value)
        
        return input_value



def get_output_from_lookup_table(file_path, input_value, temperature):
    # Load the lookup table (assuming CSV format with headers)
    lookup_table = pd.read_csv(file_path)
    
    # Extract column names
    headers = lookup_table.columns.tolist()
    
    # First column: input values
    inputs = lookup_table[headers[0]].values
    
    # Temperature steps derived from column names (e.g., "Output at 0C", "Output at 10C", ...)
    temperatures = [int(header.split(" ")[-1][:-1]) for header in headers[1:]]
    
    # Remaining columns: output data
    data = lookup_table.iloc[:, 1:].values
    
    # Interpolate/extrapolate for the input value across all temperatures
    if input_value < inputs.min() or input_value > inputs.max():
        # Extrapolation for input if outside bounds
        input_func = interp1d(inputs, data, axis=0, fill_value="extrapolate")
        outputs_at_temperatures = input_func(input_value)
    else:
        # Interpolation for input if within bounds
        input_func = interp1d(inputs, data, axis=0)
        outputs_at_temperatures = input_func(input_value)
    
    # Interpolate/extrapolate for the temperature
    if temperature < temperatures[0] or temperature > temperatures[-1]:
        # Extrapolation for temperature if outside bounds
        temp_func = interp1d(temperatures, outputs_at_temperatures, fill_value="extrapolate")
        output = temp_func(temperature)
    else:
        # Interpolation for temperature if within bounds
        temp_func = interp1d(temperatures, outputs_at_temperatures)
        output = temp_func(temperature)
    
    return output



# Example usage
file_path = "Sample_ESR_DoD.csv"  # Path to the lookup table file
input_value = 50  # Example input value
temperature = 25  # Example temperature in Celsius

# Initialize the estimator
estimator = LookupTableEstimator(file_path)

output = get_output_from_lookup_table(file_path, input_value, temperature)
print(f"Output at input {input_value} and temperature {temperature}°C: {output}")
output = estimator.output(input_value, temperature)

print(f"Output at input {input_value} and temperature {temperature}°C: {output}")

inputValue=estimator.estimate_input(output, temperature)