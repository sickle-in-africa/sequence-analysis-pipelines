---
layout: default
title: a lambda virus project
permalink: /tutorials/lambda-virus/
---

### Navigation
* [home](../../)
* [tutorials](../../tutorials/)
* [how-tos](../../how-tos/)

***

A lambda virus project
=====================================================
{: .no_toc}

Included in the base SAP library package is a reference sequence for the lambda virus, taken from the [bowite2 alignment package](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). We will use this and a tool called *simulate*, a program adapted from the simulation tool also supplied in the same bowtie2 package, to simulate reads for a fake cohort of lambda virus samples and then apply the SAP library to analyse the sequences. We will be able to check our results against the simulation's truth files, and compare several different pipelines. This tutorial will get you up and running with the library. 

In this quick tutorial we will learn how to:
* install the SAP library using the setup script
* simulate fake reads for an imaginary cohort with the `simulate` tool
* joint-call the variants from this simulated data set

[Next](1_setup.md)

### tutorial steps
1. [setting up the library](1_setup.md)
1. [simulating reads](2_simulate.md)
1. [calling variants](3_call-variants.md)
