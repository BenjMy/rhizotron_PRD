# rhizotron PRD
This repo contains the data and codes to reproduce results from a laboratory PRD experiment. 
 
## ERT and MALM
Basic import and processing using cmd line

``` {bash}
from pyPRD import processing as proc
import surveyPRD
```

- All cycle ERT analysis 
``` {bash}
python invert.py -cycle 0 1 2 3 4 5 6 7 8 9 -TL 0 -icsd 0 -reprocessed 0 -filter_seq 1 -filter_seq_rec 0 -recErr 5
```

-  ICSD analysis 
``` {bash}
python invert.py -cycle 3 4 5 6 7 8 9 -TL 0 -icsd 1 -reprocessed 0 -filter_seq 1 -filter_seq_rec 0 -recErr 5
```


## Load cells
``` {bash}
surveyPRD.load_scale_data(args.startD,args.endD)
```

## Data assimilation





# Python package for ERT and MALM data processing


### Dependencies
- Resipy
- Icsd
- pyCATHY 

### Installation and testing

* Create a conda env
- python 3.9
* Install using setup.py
!git clone the github pyPRD-dev repo
python setup.py develop or install

+ install icsd3d

for dev version:

!git clone the github icsd-dev repo
python setup.py develop or install

* Or with pip: "pip install pyPRD" [**NOT YET DEPLOYED**]


---

### Licence ###
This is free software: you can redistribute it and/or modify it under the terms of the BSD 3-clause License. A copy of this license is provided in LICENSE.md.

### References ###



