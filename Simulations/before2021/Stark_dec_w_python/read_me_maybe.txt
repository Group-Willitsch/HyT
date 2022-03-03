This is a read me file for the simulations.

Status Mid Nov 2020 (also commented in the main.py):
- You can calculate a guiding sequence and use it for the experiment
- You can calculate a decelerating sequence. It is--however--not yet perfectly matching with the sequences from the fortran code.
This could be due to:
        - A different integrator, i.e., Runge Kutta implementation missing. Currently Verlet/Euler.
        - Other interpolator for the acceleration fields
        - A correction from for the deceleration I heard from Claudio: The electrodes of the decelerator have different width.
        The were sorted in the experiment by either increasing or decreasing diameter. This apparently had to be taken into account do get good results for the sequences.
- Calculation of a decelerated molecular beam has still a lot of problems, just look at it ^^
- The code is horribly slow due due to the interpolation and integration process.
  - find_gradient_ast and find_gradient_ast_2 were first tries to find the gradient fast. I vaguely remember that I was close but that it always had problems

How to:

Guiding/deceleration sequence calculation:
- make sure you calculate only one molecule self.N = 1e0, that the molecule starts in the center etc pp.
- time steps for proper sequence should be at least self.dt = 5 ns
- Set the variables
self.sequence_calculation = True                #####If true calculate switching sequence
self.save_sequence_for_experiment = True        #####If true save the calculated sequence
- With self.phase_deg = XYZ.Y you can determine the phase angle -> the final velocity
- Search for "self.save_sequence_for_experiment" to figure out where the sequence is saved.
- In OneNote/HyT/Thomas/Simulations Experiment/Generating Dec-Sequence for experiment/Adapting the experimental code
you can find how you can use the calculated sequence in the experiment (pulseseq.py)


Calculate molecular beam with given guiding sequence:
- set self.sequence_calculation = False
- Make sure you have a guiding sequence and and read it in in:
self.read_switching_cycle_list = np.load('switching/switching_cycle_sim_10kV_5e-07_0.0deg.npy')
NOTE: TIME STEPS SIZE IS IMPORTANT HERE AND HAS TO MATCH UR CURRENT SETTINGS
- Select your molecular beam what you want to calculate, i.e.,
