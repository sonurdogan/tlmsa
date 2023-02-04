# tlmsa
Detection of SUMOylation sites that emerge through mutations in cancer. The pipeline uses the [`SUMOnet`](https://github.com/berkedilekoglu/SUMOnet) to predict possible SUMOylation sites.

## Running the Pipeline
Pipeline consist of three main part:

- Retrieving Mutation data from GDC Database and filtering respected to patients gene that has mutation resulted in lysine and getting all of the mutations of the corresponding genes of the patient since mutation near the mutated K may affect SUMOylation (R code).
- Getting wild type sequence and mapping the mutations to wild-type sequence to have mutated sequence of each protein of patient.
- Constructing 21 long subsequence (mutated K in the middle) as a input of SUMOnet.

Pipeline can be performed using bash script with a input of path and project name.

```shell
./bash_script/tlmsa.sh
```


Part 1 can be done by retrieveData.R by defining TCGA project name in the code. Once data is retrieved, Part 2 and 3 can be found as a part of tlmsa python package. 

Also, Part 2 and 3 can perform on a data other than TCGA. Detailed instraction can be found in tutorial.py 

-Pipeline workflow-


 ![workflow](https://user-images.githubusercontent.com/72014272/216789003-93ad3991-2f1f-44da-a028-f084efab50bc.jpg)



