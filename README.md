# PlastidHub<br />
An integrated platform with tools for analyzing plastome sequences

# Get Started
PlastidHub is an integrated platform with developed or collected open source tools for analyzing plastome sequences.

There are a comprehensive list of tools, from the initial sequence annotation to the final phylogenetic and comparative genomic analysis. Ultimately it will provide assistance to academic users.

Currently, PlastidHub can perform annotation, assessment of annotation quality, submission of sequence to GenBank, visualization, extraction of sequence, processing before alignment (pre-alignment), sequence alignment and processing after alignment, phylogenetic reconstruction, and comparative genomics for plastome sequences.

# User Guide
If you are new to PlastidHub, we recommend you to read through document of each tool to learn more about its functionality.

In general, wrongly changing parameter value will show an error warning alerting you that you should choose a correct value within the stated threshold range, and wrongly uploading files will show an error warning alerting you that you should upload file(s) with correct format, restricted number, and restricted size.

# Tech Stack
It is a web application platform developed with the Java,TypeScript and Perl, which is comfortable for users who prefer browser-driven applications or those who are not familiar with UNIX command-line tools.

Task-based user interface design is the most important feature of PlastidHub. The front-end web interface was developed with the vue3 (https://cn.vuejs.org/) framework and Element-Plus (https://element-plus.org/).

The back-end web framework was SpringBoot v3.3.2 (https://docs.spring.io/spring-boot/index.html) with RabbitMQ as message queue and Redis stores result status, which can manage data and parameters submitted by users, and then relay the analysis results (generated files) to the front-end for users to manage (select all, clear selection, download, and delete).

The back-end tools of PlastidHub were developed with Perl v5.26.1 programming language. The dependent third-party tools are BLAST+ v2.16.0, Circos v0.66, MAFFT v7.3, and RAxML v8.2.10.

# Browser Compatibility
PlastidHub has been tested on recently released web browsers, such as Google Chrome, Microsoft Edge, Mozilla Firefox, Opera, Apple Safari, 360 Speed Browser X, 360 Secure Browser, QQ Browser, and Huawei Browser.
