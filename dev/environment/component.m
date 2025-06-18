%{
import openmdao.api as om
import numpy as np
}%

%{
class environmentComponent(om.ExplicitComponent):
    def setup(self):

        #adding inputs:
        some_float_value = 0.0
        self.add_input(name = 'firstFloatInput', val = 0.0, desc="this is my first float variable")
        self.add_input(name = 'firstFloatMatrixInput', val = np.zeros((4,5)), shape = (4,5), desc="this is my first float matrix variable" )

        #adding output:
        self.add_output(name = 'firstFloatOutput', desc = "this is my first float variable output")
        self.add_output(name = 'firstFloatMatrixOutput', shape=(4,5), desc = "this is my first float matrix variable output")

        # Partial derivatives required for optimization
        self.declare_partials('*', '*', method='fd')


    def compute(self, inputs, outputs):
        #retrieve inputs
        firstFloatInput = inputs['firstFloatInput'][0] #[0] index, retrieve float number, otherwise it is an array
        firstFloatMatrixInput = inputs['firstFloatMatrixInput']

        firstFloatOutput = 2 * firstFloatInput
        firstFloatMatrixOutput = firstFloatOutput * firstFloatMatrixInput

        #assign outputs
        outputs['firstFloatOutput'] = firstFloatOutput
        outputs['firstFloatMatrixOutput'] = firstFloatMatrixOutput
}%

%{
#componentTest
prob = om.Problem()
#subsystem name as test
prob.model.add_subsystem('test', environmentComponent())
prob.setup()
firstFloatInput = 0.2
firstFloatMatrixInput = np.ones((4,5))
}%


&{
prob.set_val('test.firstFloatInput', firstFloatInput)
prob.set_val('test.firstFloatMatrixInput', firstFloatMatrixInput)

prob.run_model()

prob.model.list_inputs()
"""
varname                  val                 prom_name                 
-----------------------  ------------------  --------------------------
test
  firstFloatInput        [0.2]               test.firstFloatInput      
  firstFloatMatrixInput  |4.47213595|        test.firstFloatMatrixInput
"""

prob.model.list_outputs()
"""
varname                   val                   prom_name                  
------------------------  --------------------  ---------------------------
test
  firstFloatOutput        [0.4]                 test.firstFloatOutput      
  firstFloatMatrixOutput  |1.78885438|          test.firstFloatMatrixOutput
"""

print(prob.get_val('test.firstFloatMatrixOutput')) # to get matrix output

"""
[[0.4 0.4 0.4 0.4 0.4]
 [0.4 0.4 0.4 0.4 0.4]
 [0.4 0.4 0.4 0.4 0.4]
 [0.4 0.4 0.4 0.4 0.4]]


"""
 
}%

% Plan: rewrite everything using the optimization toolboxes in MATLAB.


