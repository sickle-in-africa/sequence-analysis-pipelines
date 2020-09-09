---
layout: default
title: clean downstream data
parent: how-to's
---

How do I quickly clean data downstream of read files?
======================================================

Sometimes we make a mistake, and the pipeline prematurely terminates or produces tainted or even non-sensical data. If we have just encountered such a scenario, and we want to clean all the pipeline data *downstream of the raw read files* we can use the clean-data pipeline for whole genome sequence data.

## tools required

just the SAP library basic installation

## steps

1. make a note of the cohort ID for the faulty run, for example `c_2`. Create a json file in the inputs folder, `<sap-path>/inputs/wgs` that contains the cohort id:
```
{
  "cohort_id": "c_1"
}
```
and save it as:
```
	<sap-path>/inputs/wgs/clean-data.c_1.json
```
where in the above, `c_1` should be changed to reflect the id of the cohort whose downstream data you wish to flush

1. Change to the root directory, `<sap-path>`, and run the clean-data pipeline using the sap script:
```
./sap.sh wgs clean-data c_1
```

## comments
1. be very careful with this command, because it is *irreversible*. Once deleted, the downstream data cannot be retreived.
1. the raw reads files themselves are not removed, however the trimmed read files are.  