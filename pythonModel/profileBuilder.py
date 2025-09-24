# -*- coding: utf-8 -*-
"""
Created on Sun Apr 13 10:27:57 2025

@author: iitm9
"""

import pandas as pd
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from typing import List, Union, Optional

class BatteryProfileSimulator:
    def __init__(self, profile_directory: str = '.', default_capacity: float = 10.0):
        """
        Initialize the battery profile simulator.
        
        Args:
            profile_directory: Directory containing the CSV profile files
            default_capacity: Default battery capacity in Amp-hours (Ah)
        """
        self.profile_directory = profile_directory
        self.profiles = {}
        self.default_capacity = default_capacity
        self.load_profiles()
    
    def load_profiles(self):
        """Load all CSV profiles from the profile directory"""
        profile_files = glob.glob(os.path.join(self.profile_directory, '*.csv'))
        
        for file_path in profile_files:
            profile_name = os.path.basename(file_path).replace('.csv', '')
            try:
                df = pd.read_csv(file_path)
                if all(col in df.columns for col in ['Time', 'Current', 'Temperature']):
                    self.profiles[profile_name] = df
                    print(f"Loaded profile: {profile_name}")
                else:
                    print(f"Warning: {profile_name} is missing required columns")
            except Exception as e:
                print(f"Error loading {profile_name}: {e}")
        
        print(f"Loaded {len(self.profiles)} profiles")
    
    def get_available_profiles(self):
        """Return a list of available profile names"""
        return list(self.profiles.keys())
    
    def scale_profile_current(self, profile: pd.DataFrame, target_capacity: float) -> pd.DataFrame:
        """
        Scale the current values in a profile to match a target battery capacity.
        
        Args:
            profile: DataFrame containing the profile data
            target_capacity: Target battery capacity in Amp-hours
            
        Returns:
            DataFrame with scaled current values
        """
        if target_capacity <= 0:
            raise ValueError("Target capacity must be positive")
            
        # Create a copy to avoid modifying the original
        scaled_profile = profile.copy()
        
        # Scale current based on the ratio of target capacity to default capacity
        scaling_factor = target_capacity / self.default_capacity
        scaled_profile['Current'] = profile['Current'] * scaling_factor
        
        return scaled_profile
    
    def create_composite_profile(self, profile_names: List[str], cycles: int = 1, 
                                capacity: Optional[float] = None) -> pd.DataFrame:
        """
        Create a composite profile by concatenating selected profiles.
        
        Args:
            profile_names: List of profile names to include
            cycles: Number of times to repeat the sequence
            capacity: Target battery capacity in Amp-hours (if None, use defaults)
            
        Returns:
            DataFrame with the composite profile
        """
        if not profile_names:
            raise ValueError("No profiles selected")
        
        # Verify all profiles exist
        missing = [name for name in profile_names if name not in self.profiles]
        if missing:
            raise ValueError(f"Profiles not found: {', '.join(missing)}")
        
        # Create cycle of selected profiles
        cycle_dfs = []
        
        for _ in range(cycles):
            for name in profile_names:
                df = self.profiles[name].copy()
                
                # Scale current if a capacity is specified
                if capacity is not None and capacity != self.default_capacity:
                    df = self.scale_profile_current(df, capacity)
                
                # If not the first profile in the sequence, adjust the time
                if cycle_dfs:
                    last_time = cycle_dfs[-1]['Time'].iloc[-1]
                    df['Time'] = df['Time'] + last_time
                
                cycle_dfs.append(df)
        
        # Concatenate all dataframes
        result = pd.concat(cycle_dfs, ignore_index=True)
        return result
    
    def calculate_energy_metrics(self, profile: pd.DataFrame) -> dict:
        """
        Calculate energy metrics for a profile.
        
        Args:
            profile: DataFrame containing the profile data
            
        Returns:
            Dictionary with energy metrics
        """
        metrics = {}
        
        # Split into charge and discharge
        discharge = profile[profile['Current'] > 0].copy()
        charge = profile[profile['Current'] < 0].copy()
        
        # Calculate amp-hours and energy
        discharge_ah = 0
        charge_ah = 0
        
        # Nominal voltage (assumed)
        nominal_voltage = 3.7  # Li-ion typical
        
        # Calculate discharge amp-hours
        if len(discharge) > 1:
            for i in range(1, len(discharge)):
                current_avg = (discharge['Current'].iloc[i-1] + discharge['Current'].iloc[i]) / 2
                time_diff = (discharge['Time'].iloc[i] - discharge['Time'].iloc[i-1]) / 3600  # convert to hours
                discharge_ah += current_avg * time_diff
        
        # Calculate charge amp-hours
        if len(charge) > 1:
            for i in range(1, len(charge)):
                current_avg = abs((charge['Current'].iloc[i-1] + charge['Current'].iloc[i]) / 2)
                time_diff = (charge['Time'].iloc[i] - charge['Time'].iloc[i-1]) / 3600  # convert to hours
                charge_ah += current_avg * time_diff
        
        # Calculate watt-hours
        discharge_wh = discharge_ah * nominal_voltage
        charge_wh = charge_ah * nominal_voltage
        
        # Save metrics
        metrics['discharge_ah'] = discharge_ah
        metrics['charge_ah'] = charge_ah
        metrics['discharge_wh'] = discharge_wh
        metrics['charge_wh'] = charge_wh
        metrics['duration_seconds'] = profile['Time'].max() - profile['Time'].min()
        metrics['duration_hours'] = metrics['duration_seconds'] / 3600
        
        return metrics
    
    def simulate(self, profile_names: List[str], cycles: int = 1, 
                capacity: Optional[float] = None, output_file: Optional[str] = None, 
                visualize: bool = False) -> pd.DataFrame:
        """
        Run a simulation with the selected profiles for a specified number of cycles.
        
        Args:
            profile_names: List of profile names to include
            cycles: Number of times to repeat the sequence
            capacity: Target battery capacity in Amp-hours (if None, use defaults)
            output_file: If provided, save the composite profile to this CSV file
            visualize: Whether to visualize the profile
            
        Returns:
            DataFrame with the composite profile
        """
        # Use default capacity if none specified
        if capacity is None:
            capacity = self.default_capacity
            
        # Create the composite profile
        composite = self.create_composite_profile(profile_names, cycles, capacity)
        
        # Calculate energy metrics
        metrics = self.calculate_energy_metrics(composite)
        
        # Print summary
        print(f"Created composite profile with {len(composite)} time points")
        print(f"Battery capacity: {capacity} Ah")
        print(f"Duration: {metrics['duration_hours']:.2f} hours")
        print(f"Total discharge: {metrics['discharge_ah']:.2f} Ah ({metrics['discharge_wh']:.2f} Wh)")
        print(f"Total charge: {metrics['charge_ah']:.2f} Ah ({metrics['charge_wh']:.2f} Wh)")
        
        # Verify that discharge doesn't exceed battery capacity
        if metrics['discharge_ah'] > capacity:
            print(f"WARNING: Total discharge ({metrics['discharge_ah']:.2f} Ah) exceeds "
                  f"battery capacity ({capacity} Ah)")
        
        # Save to file if requested
        if output_file:
            composite.to_csv(output_file, index=False)
            print(f"Saved composite profile to {output_file}")
        
        # Visualize if requested
        if visualize:
            self.visualize_profile(composite, capacity)
        
        return composite
    
    def visualize_profile(self, profile: pd.DataFrame, capacity: float = None):
        """
        Visualize the given profile with current, temperature, and SOC plots.
        
        Args:
            profile: DataFrame containing the profile data
            capacity: Battery capacity for reference line
        """
        if capacity is None:
            capacity = self.default_capacity
            
        # Calculate SOC estimation
        soc = self.estimate_soc(profile, capacity)
        
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
        
        # Current plot
        ax1.plot(profile['Time'], profile['Current'], 'b-')
        ax1.axhline(y=0, color='k', linestyle='-', alpha=0.3)
        ax1.fill_between(profile['Time'], profile['Current'], 0, 
                         where=(profile['Current'] < 0), color='g', alpha=0.3, 
                         label='Charging')
        ax1.fill_between(profile['Time'], profile['Current'], 0, 
                         where=(profile['Current'] > 0), color='r', alpha=0.3, 
                         label='Discharging')
        ax1.set_ylabel('Current (A)')
        ax1.set_title('Battery Current Profile')
        ax1.legend()
        ax1.grid(True)
        
        # Temperature plot
        ax2.plot(profile['Time'], profile['Temperature'], 'r-')
        ax2.set_ylabel('Temperature (°C)')
        ax2.set_title('Ambient Temperature Profile')
        ax2.grid(True)
        
        # SOC plot
        ax3.plot(profile['Time'], soc, 'g-')
        ax3.set_xlabel('Time (seconds)')
        ax3.set_ylabel('Estimated SOC (%)')
        ax3.set_title('Estimated State of Charge')
        ax3.set_ylim(0, 100)
        ax3.grid(True)
        
        plt.tight_layout()
        plt.show()
    
    def estimate_soc(self, profile: pd.DataFrame, capacity: float = None) -> pd.Series:
        """
        Estimate SOC based on current integration (simplified).
        
        Args:
            profile: DataFrame containing the profile data
            capacity: Battery capacity in Ah
            
        Returns:
            Series with estimated SOC values
        """
        if capacity is None:
            capacity = self.default_capacity
            
        # Initialize SOC at 50% (arbitrary starting point)
        initial_soc = 50
        soc = np.zeros(len(profile))
        soc[0] = initial_soc
        
        # Coulomb counting (simple integration)
        for i in range(1, len(profile)):
            # Current average in this time step
            current_avg = (profile['Current'].iloc[i-1] + profile['Current'].iloc[i]) / 2
            
            # Time difference in hours
            time_diff = (profile['Time'].iloc[i] - profile['Time'].iloc[i-1]) / 3600
            
            # Change in SOC (negative current means charging)
            delta_soc = -100 * (current_avg * time_diff) / capacity
            
            # Update SOC
            soc[i] = soc[i-1] + delta_soc
            
        # Clip SOC to valid range [0, 100]
        soc = np.clip(soc, 0, 100)
        
        return pd.Series(soc)


# Example usage
if __name__ == "__main__":
    # Initialize the simulator
    simulator = BatteryProfileSimulator('usageProfile')
    #simulator.load_profiles()
    
    # Get available profiles
    profiles = simulator.get_available_profiles()
    print(f"Available profiles: {profiles}")
    

    
    # Example 3: Create a custom test sequence
    test_sequence = [
        'EVdeepDischarge',  # Deep discharge
        'EVfastCharge',  # Fast charge
        'EVhighwayCommute',  # Highway trip
        'EVpartialCharge',  # Partial charge
        'EVregularCommute'
    ]
    
    profile3 = simulator.simulate(
        profile_names=test_sequence,
        cycles=10,
        capacity=65.0,
        output_file='test_sequence.csv',
        visualize=True
    )