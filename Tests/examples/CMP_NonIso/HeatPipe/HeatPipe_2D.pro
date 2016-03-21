<?xml version="1.0"?>
<ogs6>
<coupling>
    <P algorithm="Serial" convergence="FemFunctionConvergenceCheck" max_itr="1" epsilon="1e-4">
        <out>MEAN_PRESSURE</out>
        <out>TOTAL_MASS_DENSITY</out>
		<out>TEMPERATURE</out>
        <problems>
            <M name="CMP_NonIso_TotalDensityForm" type="CMP_NonIso_TotalDensityForm">
                <out>MEAN_PRESSURE</out>
				<out>TOTAL_MASS_DENSITY</out>
				<out>TEMPERATURE</out>
            </M>
        </problems>
    </P>
</coupling>
</ogs6> 

