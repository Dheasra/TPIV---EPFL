import math
import numpy as np
from numpy import pi
# importing Qiskit
from qiskit import QuantumCircuit, execute, Aer, IBMQ, ClassicalRegister, QuantumRegister
from qiskit.providers.ibmq import least_busy
from qiskit.providers.aer import QasmSimulator
from qiskit.tools.monitor import job_monitor
from qiskit.visualization import plot_histogram, plot_bloch_multivector
from qiskit.providers.aer.noise import NoiseModel
from qiskit.providers.aer.noise.errors import pauli_error, depolarizing_error
from qiskit.compiler import transpile 
# from qiskit.compiler.assemble import assemble
# from qiskit.assembler.disassemble import disassemble
# %config InlineBackend.figure_format = 'svg' # Makes the images look nice
import matplotlib.pyplot as plt
import matplotlib

def get_noise(p, qubits=1): #from: https://github.com/MIGUEL-LO/Qiskit_depolarazation_channel/blob/master/depo_channel_using_depolarizing_error.ipynb
    
    # This creates the depolarizing error channel,
    # epsilon(P) = (1-P)rho + (P/3)(XrhoX + YrhoY + ZrhoZ).
    depo_err_chan = depolarizing_error((4*p)/3, qubits)

    # Creating the noise model to be used during execution.
    noise_model = NoiseModel()

    noise_model.add_all_qubit_quantum_error(depo_err_chan, "measure") # measurement error is applied to measurements

    return noise_model, depo_err_chan


qreg_q = QuantumRegister(3, 'q')
creg_c = ClassicalRegister(3, 'c')
circuit = QuantumCircuit(qreg_q, creg_c)

circuit.u2(0, pi, qreg_q[0])
circuit.u2(0, pi, qreg_q[1])
circuit.u2(0, pi, qreg_q[2])
circuit.barrier(qreg_q[0])
circuit.barrier(qreg_q[1])
circuit.barrier(qreg_q[2])
circuit.u(-pi, pi/2, 3*pi/2, qreg_q[0])
circuit.u(pi, 0, pi, qreg_q[1])
circuit.u(pi, 0, pi, qreg_q[2])
circuit.cx(qreg_q[0], qreg_q[1])
circuit.id(qreg_q[2])
circuit.cx(qreg_q[1], qreg_q[0])
circuit.id(qreg_q[2])
circuit.cx(qreg_q[0], qreg_q[1])
circuit.id(qreg_q[2])
circuit.id(qreg_q[0])
circuit.cx(qreg_q[2], qreg_q[1])
circuit.id(qreg_q[0])
circuit.u1(-pi/4, qreg_q[1])
circuit.id(qreg_q[2])
circuit.cx(qreg_q[0], qreg_q[1])
circuit.id(qreg_q[2])
circuit.id(qreg_q[0])
circuit.u1(pi/4, qreg_q[1])
circuit.id(qreg_q[2])
circuit.id(qreg_q[0])
circuit.cx(qreg_q[2], qreg_q[1])
circuit.id(qreg_q[0])
circuit.u1(-pi/4, qreg_q[1])
circuit.u1(pi/4, qreg_q[2])
circuit.cx(qreg_q[1], qreg_q[0])
circuit.id(qreg_q[2])
circuit.cx(qreg_q[0], qreg_q[1])
circuit.id(qreg_q[2])
circuit.u2(-pi, 13.4, qreg_q[0])
circuit.cx(qreg_q[1], qreg_q[2])
circuit.id(qreg_q[0])
circuit.u1(pi/4, qreg_q[1])
circuit.u1(-pi/4, qreg_q[2])
circuit.id(qreg_q[0])
circuit.cx(qreg_q[1], qreg_q[2])
circuit.id(qreg_q[0])
circuit.u2(-pi, 2*pi, qreg_q[1])
circuit.u2(-pi, 2*pi, qreg_q[2])
circuit.id(qreg_q[0])
circuit.cx(qreg_q[1], qreg_q[2])
circuit.id(qreg_q[0])
circuit.cx(qreg_q[2], qreg_q[1])
circuit.id(qreg_q[0])
circuit.cx(qreg_q[1], qreg_q[2])
circuit.cx(qreg_q[1], qreg_q[0])
circuit.id(qreg_q[2])
circuit.u1(-pi/4, qreg_q[0])
circuit.cx(qreg_q[1], qreg_q[2])
circuit.id(qreg_q[0])
circuit.cx(qreg_q[2], qreg_q[1])
circuit.id(qreg_q[0])
circuit.cx(qreg_q[1], qreg_q[2])
circuit.cx(qreg_q[1], qreg_q[0])
circuit.id(qreg_q[2])
circuit.u1(pi/4, qreg_q[0])
circuit.cx(qreg_q[1], qreg_q[2])
circuit.id(qreg_q[0])
circuit.cx(qreg_q[2], qreg_q[1])
circuit.id(qreg_q[0])
circuit.cx(qreg_q[1], qreg_q[2])
circuit.cx(qreg_q[1], qreg_q[0])
circuit.id(qreg_q[2])
circuit.u1(-pi/4, qreg_q[0])
circuit.u1(pi/4, qreg_q[1])
circuit.id(qreg_q[2])
circuit.id(qreg_q[0])
circuit.cx(qreg_q[1], qreg_q[2])
circuit.id(qreg_q[0])
circuit.cx(qreg_q[2], qreg_q[1])
circuit.id(qreg_q[0])
circuit.cx(qreg_q[1], qreg_q[2])
circuit.cx(qreg_q[1], qreg_q[0])
circuit.id(qreg_q[2])
circuit.u2(pi, 13*pi/4, qreg_q[0])
circuit.cx(qreg_q[1], qreg_q[2])
circuit.id(qreg_q[0])
circuit.u1(pi/4, qreg_q[1])
circuit.u1(-pi/4, qreg_q[2])
circuit.id(qreg_q[0])
circuit.cx(qreg_q[1], qreg_q[2])
circuit.id(qreg_q[0])
circuit.u2(pi, pi, qreg_q[1])
circuit.u2(pi, pi, qreg_q[2])
circuit.barrier(qreg_q[0])
circuit.barrier(qreg_q[1])
circuit.barrier(qreg_q[2])
circuit.measure(qreg_q[0], creg_c[0])
circuit.measure(qreg_q[1], creg_c[1])
circuit.measure(qreg_q[2], creg_c[2])


#visualization of the circuit
# print(grvr) #in the terminal
circuit.draw('mpl')

#control panel
ns = 1 #1 = sim w/ noise, 2 = real device, else = sim
n_model = 0 #0 = model from backend, 1= depolarizing error
P = 0.2 #depolarizing probability

#=== Results ===
#Getting the noise model
shots = 2048 #nbr of runs of the circuit

if ns == 1: #sim with noise
    provider         = IBMQ.load_account()

    backend          = provider.get_backend('ibmq_valencia')
    if n_model == 0:
        noise_model  = NoiseModel.from_backend(backend)
    else:
        noise_model, depo_err_chan = get_noise(P, 1)
    
    simulator        = Aer.get_backend('qasm_simulator')
    results = execute(circuit, backend=simulator, shots=shots,noise_model = noise_model).result()
elif ns == 2: #real device
    provider         = IBMQ.load_account()
    backend         = provider.get_backend('ibmq_valencia')
    job = execute(circuit, backend = backend, shots = 1024, optimization_level= 3)
    job_monitor(job, interval= 2)
    results = job.result()

else: #sim
    backend = Aer.get_backend('qasm_simulator') #selection of the device on which to execute the circuit
    results = execute(circuit, backend = backend, shots = shots).result()


# qobj = assemble(circuit, backend, shots)
# comp_circ = disassemble(qobj)
# comp_circ = comp_circ[0][0]
# comp_circ.draw('mpl')

#Getting the results
answer = results.get_counts()

#ploting
plot_histogram(answer)
plt.show()