#! /usr/bin/env python
"""
Script for generation of input files for Aspirin and Naproxen binding to HSA by ITC.
"""

from itctools.itctools import ureg, Quantity
from itctools.procedures import ITCProtocol, ITCExperimentSet, ITCExperiment, ITCHeuristicExperiment
from itctools.materials import Solvent, Compound, SimpleSolution
from itctools.labware import Labware, PipettingLocation
import pprint


# The sample cell volume in microliters
cell_volume = 202.8

# Define solvents.
water = Solvent('water', density=0.9970479 * ureg.gram / ureg.milliliter)
buffer = Solvent('buffer', density=1.014 * ureg.gram / ureg.milliliter) # TODO is our density the same as the HOST-GUEST buffer?

# Define compounds.
hsa = Compound('HumanSerumAlbumin', molecular_weight=65000 * (ureg.gram / ureg.mole), purity=.95)
naproxen_sodium = Compound('Naproxen Sodium', molecular_weight=252.24 * (ureg.gram / ureg.mole), purity=1.)
# TODO : update aspirin molecular weight - MI
aspirin = Compound('AcetylsalicylicAcid', molecular_weight=180.15742 * (ureg.gram / ureg.mole), purity=.99)


# Ka (association constants) TODO Add this to the compound properties? (maybe a dict with protein as key)
aspirin_ka = 2.76 / ureg.micromolar  # http://omicsonline.org/2157-7544/2157-7544-2-107.pdf (first site estimate at 298K)
naproxen_ka = 1.73e7 / ( ureg.molar)  # http://pubs.acs.org/doi/pdf/10.1021/jp062734p with (BSA, with 200 uM NaCl)

# Define troughs on the instrument
water_trough = Labware(RackLabel='Water', RackType='Trough 100ml')
buffer_trough = Labware(RackLabel='Buffer', RackType='Trough 100ml')

# Define source labware.
source_plate = Labware(RackLabel='SourcePlate', RackType='5x3 Vial Holder')

# Define source solutions in the vial holder
#TODO : Define solutions once prepared with the Quantos. Enter compound mass and solvent mass values. - MI
#I entered compound mass and solvent mass that matches  40 uM
hsa_solution = SimpleSolution(compound=hsa, compound_mass=9.40 * ureg.milligram, solvent=buffer, solvent_mass=3.4245 * ureg.gram, location=PipettingLocation(
    source_plate.RackLabel,
    source_plate.RackType,
    1))

aspirin_solution = SimpleSolution(compound=aspirin, compound_mass= 20.235 * ureg.milligram, solvent=buffer, solvent_mass= 10.0973 * ureg.gram, location=PipettingLocation(
    source_plate.RackLabel,
    source_plate.RackType,
    2))

naproxen_sodium_solution = SimpleSolution(compound=naproxen_sodium, compound_mass=25.275 * ureg.milligram, solvent=buffer, solvent_mass=10.0845 * ureg.gram, location=PipettingLocation(
    source_plate.RackLabel,
    source_plate.RackType,
    3))


drugs = [aspirin, naproxen_sodium]
drug_solutions = [aspirin_solution, naproxen_sodium_solution]
drug_kas = [aspirin_ka, naproxen_ka]



# Define ITC protocol.

# Protocol for 'control' titrations (water-water, buffer-buffer,
# titrations into buffer, etc.)

control_protocol = ITCProtocol(
    'control_protocol',
    sample_prep_method='Chodera Load Cell Without Cleaning Cell After.setup',
    itc_method='ChoderaWaterWater5.inj',
    analysis_method='Control',
    experimental_conditions=dict(target_temperature=25, equilibration_time=60, stir_rate=1000, reference_power=5),
    injections=[
        dict(volume_inj=0.1, duration_inj=0.4, spacing=60, filter_period=0.5)] +
        5 * [dict(volume_inj=1.5, duration_inj=6, spacing=120, filter_period=0.5)],
    )


blank_protocol = ITCProtocol(
    '1:1 binding protocol',
    sample_prep_method='Chodera Load Cell Without Cleaning Cell After.setup',
    itc_method='ChoderaHSA20.inj',
    analysis_method='Onesite',  # Protocol for 1:1 binding analysis
    experimental_conditions=dict(target_temperature=25, equilibration_time=300, stir_rate=1000, reference_power=5),
    injections=[
        dict(volume_inj=0.1, duration_inj=0.4, spacing=60, filter_period=0.5)] +
        20 * [dict(volume_inj=1.5, duration_inj=6, spacing=120, filter_period=0.5)],
    )


binding_protocol = ITCProtocol(
    '1:1 binding protocol',
    sample_prep_method='Plates Standard.setup', # includes cleaning in the end
    itc_method='ChoderaHSA20.inj',
    analysis_method='Onesite',
    experimental_conditions=dict(target_temperature=25, equilibration_time=300, stir_rate=1000, reference_power=5),
    injections=[
        dict(volume_inj=0.1, duration_inj=0.4, spacing=60, filter_period=0.5)] +
        20 * [dict(volume_inj=1.5, duration_inj=6, spacing=120, filter_period=0.5)],
    )


cleaning_protocol = ITCProtocol(
    'cleaning protocol',
    sample_prep_method='Plates Clean.setup',
    itc_method='ChoderaWaterWater5.inj',
    analysis_method='Onesite',
    experimental_conditions=dict(target_temperature=25, equilibration_time=60, stir_rate=1000, reference_power=5),
    injections=5 * [
        dict(volume_inj=7.5, duration_inj=15, spacing=150, filter_period=5)],
    )

# Define ITC Experiment.

# use specified protocol by default
itc_experiment_set = ITCExperimentSet(name='Human Serum Albumin experiments')
# Add available plates for experiments.
itc_experiment_set.addDestinationPlate(
    Labware(
        RackLabel='DestinationPlate',
        RackType='ITC Plate'))
itc_experiment_set.addDestinationPlate(
    Labware(
        RackLabel='DestinationPlate2',
        RackType='ITC Plate'))

nreplicates = 1  # number of replicates of each experiment


# #  Initial cleaning is skipped.
# name = 'initial cleaning water titration'
# itc_experiment_set.addExperiment(
#     ITCExperiment(
#         name=name,
#         syringe_source=water_trough,
#         cell_source=water_trough,
#         protocol=cleaning_protocol,
#         cell_volume=cell_volume))

# Add water control titrations.
for replicate in range(1):
    name = 'water into water %d' % (replicate + 1)
    itc_experiment_set.addExperiment(
        ITCExperiment(
            name=name,
            syringe_source=water_trough,
            cell_source=water_trough,
            protocol=control_protocol,
            cell_volume=cell_volume))

# Add buffer control titrations.
for replicate in range(1):
    name = 'buffer into buffer %d' % (replicate + 1)
    itc_experiment_set.addExperiment(
        ITCExperiment(
            name=name,
            syringe_source=buffer_trough,
            cell_source=buffer_trough,
            protocol=blank_protocol,
            cell_volume=cell_volume))



# drugs/HSA
# scale cell concentration to fix necessary syringe concentrations

cell_scaling = 1.
for drug, drug_solution, drug_ka in zip(drugs, drug_solutions, drug_kas):

    # We need to store the experiments before adding them to the set
    drug_protein_experiments = list()
    drug_buffer_experiments = list()

    # Scaling factors per replicate
    factors = list()


    # Define drug to protein experiments.
    for replicate in range(1):
        name = '%s into HSA %d' % (drug.name, replicate + 1 )
        experiment = ITCHeuristicExperiment(
            name=name,
            syringe_source=drug_solution,
            cell_source=hsa_solution,
            protocol=binding_protocol,
            cell_concentration=0.040 *
            ureg.millimolar *
            cell_scaling,
            buffer_source=buffer_trough,
            cell_volume=cell_volume)
        # optimize the syringe_concentration using heuristic equations and known binding constants
        # TODO extract m, v and V0 from protocol somehow?

        ## I am not sure what this part is doing?

        # Warning, you're possibly not getting the setup you want. Consider not using the Heuristic Experiment
        experiment.heuristic_syringe(drug_ka, 10, strict=False)
        # rescale if syringe > stock. Store factor.
        factors.append(experiment.rescale())
        drug_protein_experiments.append(experiment)


    # Define drug into buffer experiments.
    for replicate in range(1):
        name = '%s into buffer  %d' % (drug.name, replicate + 1)
        experiment = ITCHeuristicExperiment(
            name=name,
            syringe_source=drug_solution,
            cell_source=buffer_trough,
            protocol=blank_protocol,
            buffer_source=buffer_trough,
            cell_volume=cell_volume)
        # rescale to match drug to protein experiment concentrations.
        experiment.rescale(tfactor=factors[replicate])
        drug_buffer_experiments.append(experiment)



    # TODO, since we are changing drugs, we'd have to wash the syringe.

    #This part determines the order of experiments

    # Add drug_to_buffer experiment(s) to set
    for drug_buffer_experiment in drug_buffer_experiments:
        itc_experiment_set.addExperiment(drug_buffer_experiment)
        # pprint.pprint(drug_buffer_experiment.__dict__)

    # Add drug to protein experiment(s) to set
    for drug_protein_experiment in drug_protein_experiments:
        itc_experiment_set.addExperiment(drug_protein_experiment)
        # pprint.pprint(drug_protein_experiment.__dict__)


# buffer into hsa
for replicate in range(1):
    name = 'buffer into HSA %d' % (replicate + 1)
    itc_experiment_set.addExperiment(
        ITCExperiment(
            name=name,
            syringe_source=buffer_trough,
            cell_source=hsa_solution,
            protocol=blank_protocol,
            cell_concentration=0.040 * ureg.millimolar,
            buffer_source=buffer_trough,
            cell_volume=cell_volume))


# Add final cleaning experiment.
name = 'final cleaning water titration'
itc_experiment_set.addExperiment(
    ITCExperiment(
        name=name,
        syringe_source=water_trough,
        cell_source=water_trough,
        protocol=cleaning_protocol,
        cell_volume=cell_volume))


# Add water control titrations.
nfinal = 1
for replicate in range(nfinal):
    name = 'final water into water test %d' % (replicate + 1)
    itc_experiment_set.addExperiment(
        ITCExperiment(
            name=name,
            syringe_source=water_trough,
            cell_source=water_trough,
            protocol=control_protocol,
            cell_volume=cell_volume))

# Check that the experiment can be carried out using available solutions
# and plates.

itc_experiment_set.validate(print_volumes=True, omit_zeroes=True)

# For convenience, concentrations
for drug_solution in drug_solutions:
    print("%s %.4f mM" % (drug_solution.name, drug_solution.concentration / ureg.millimolar ))
    print("HSA", hsa_solution.concentration.to(ureg.millimolar))


# Write Tecan EVO pipetting operations.
worklist_filename = 'hsa_20150714_MI.gwl'
itc_experiment_set.writeTecanWorklist(worklist_filename)

# Write Auto iTC-200 experiment spreadsheet.
excel_filename = 'hsa_20150714_MI.xlsx'
itc_experiment_set.writeAutoITCExcel(excel_filename)
