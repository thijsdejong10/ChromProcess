def valid_deconvolution(analysis):
    """
    Take an analysis file and verify in the deconvolution paramaters are accurate.

    """
    valid = True
    for reg in analysis.deconvolve_regions:
        sigma = analysis.deconvolve_regions[reg]["standard_sigma"]
        region_boundaries = analysis.deconvolve_regions[reg]["region_boundaries"]
        amplitude = analysis.deconvolve_regions[reg]["standard_amplitude"]
        no_peaks = analysis.deconvolve_regions[reg]["number_of_peaks"]
        if not sigma[0] <= sigma[1] <= sigma[2]:
            valid = False
            print(f"Error in {reg}: boundaries for sigma incorrect")
        if not region_boundaries[0] < region_boundaries[1]:
            valid = False
            print(f"Error in {reg}: boundaries for region incorrect")
        if not amplitude[0] <= amplitude[1] <= amplitude[2]:
            valid = False
            print(f"Error in {reg}: boundaries for amplitude incorrect")
        if not no_peaks == len(analysis.deconvolve_regions[reg]["peaks"]):
            valid = False
            print(f"Error in {reg}: boundaries for number_of_peaks incorrect")
        for peak_name, peak_dict in analysis.deconvolve_regions[reg]["peaks"].items():
            for boundaries in peak_dict:
                if (
                    not peak_dict[boundaries][0]
                    <= peak_dict[boundaries][1]
                    <= peak_dict[boundaries][2]
                ):
                    valid = False
                    print(
                        f"Error in {reg}: boundaries for {boundaries} of peak {peak_name} are incorrect."
                    )
    return valid
