# phenom-eos

Code for reproducing phenomenological models of candidate EoSs, like piecewise polytropes.

---

### Scripts

###### Piecewise polytrope (pwp)

* geteos_pwp (single EoS)
* makepwp inputparams.in (batch version)

###### Spectral EoS (spec)

* geteos_spec (single EoS)
* makespec inputparams.in (batch version)

###### Strange quark matter EoS (sqm)

* geteos_sqm (single EoS)
* makesqm inputparams.in (batch version)

###### Constant sound speed hadron-quark hybrid EoS (css)

* geteos_css (single EoS)
* makecss inputparams.in (batch version)

---

### Tools

###### EoS plot

* ploteos eos.csv

###### Affix SLy crust

* addcrust eos.csv
* addcrusts ./in/dir/ ./out/dir (batch version)

---

### References

###### Read+ PRD 79 (2009)

* Main reference for algorithm in get_pwp
* Source for parameters in ref_eos-pwp.in

###### Lindblom PRD 82 (2010)

* Main reference for algorithm in get_spec
* Source for parameters in ref_eos-spec.in

###### Han+Steiner arXiv:1810.10967

* Main reference for algorithm in get_sqm and get_css
* Source for parameters in ref_eos-sqm.in

###### Paschalidis+ PRD 97 (2018)

* Source for parameters in ref_eos-css.in and ref_eos-pwp4.in

###### Prakash+Cooke PRD 52 (1995) and Lattimer+Prakash ApJ 550 (2001)

* Secondary references for units in get_sqm

