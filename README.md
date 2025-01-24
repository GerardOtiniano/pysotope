# pysotope (v.0.1.0).

Pysotope is an open-source package meant to processes raw data measured from the OSIBL GC-IRMS. Corrections are automatically calculated but the user must verify and confirm their application.

Note: pysotope was tested using jupyter lab and jupyter notebook. Compatibility with alternative python IDEs may require custom configuration.

## Corrections

- **Drift** (weighted means least squared)
- **Linearity** (log transformed)
- **Methanol** derivitization - user is asked for methanol values
- **PAME** - calculates the isotopic composition of PAME, if run in the sequence (set argument pame=True)
- **VSMOW** - calculated using C18, C20, and C28 standards, tested on C24 standard.

## Installation

```base
pip install pysotope
```

## Arguments

- pame (default False) - if True, the package will automatically identify and correct PAME values in the sequence.
- user_linearity_conditions (default False) - if True, during linearity correction, user will be asked for a cutoff value under which samples with peak areas lower than the cutoff will be excluded from the analysis.
- alt_stds (default False) - if True, the user can modify the default standard value. These changes will be saved for the user, but will reset upon package update.

e.g.,

```bash
import pysotope
pysotope.iso_process(pame=True, user_linearity_conditions = False, alt_stds = True)
```

## Contributing

Contributions to pysotope are welcome! If you encounter any issues or have suggestions for improvements, please open an issue or submit a pull request on the GitHub repository.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions or inquiries, please contact:

- Author: Dr. Gerard Otiniano & Dr. Elizabeth Thomas
- Email: gerardot@buffalo.edu
