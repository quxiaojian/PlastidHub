# PlastidHub v1.0<br />
An integrated analysis platform for plastid phylogenomics and comparative genomics<br />
Copyright (C) 2025 Xiao-Jian Qu<br />
<br />

# Website<br />
https://www.plastidhub.cn<br />
<br />

# Contact<br />
quxiaojian@sdnu.edu.cn<br />
<br />

# Citation<br />
Na-Na Zhang, Gregory W. Stull, Xue-Jie Zhang, Shou-Jin Fan, Ting-Shuang Yi, Xiao-Jian Qu. PlastidHub: an integrated analysis platform for plastid phylogenomics and comparative genomics. Plant Diversity, 2025, https://doi.org/10.1016/j.pld.2025.05.005.<br />
<br />

# Get Started<br />
PlastidHub is an integrated platform with developed or collected open source tools for analyzing plastome sequences.<br />

There are a comprehensive list of tools, from the initial sequence annotation to the final phylogenetic and comparative genomic analysis. Ultimately it will provide assistance to academic users.<br />

Currently, PlastidHub can perform annotation, assessment of annotation quality, submission of sequence to GenBank, visualization, extraction of sequence, processing before alignment (pre-alignment), sequence alignment and processing after alignment, phylogenetic reconstruction, and comparative genomics for plastome sequences.<br />
<br />

# User Guide<br />
If you are new to PlastidHub, we recommend you to read through document of each tool to learn more about its functionality.<br />

In general, wrongly changing parameter value will show an error warning alerting you that you should choose a correct value within the stated threshold range, and wrongly uploading files will show an error warning alerting you that you should upload file(s) with correct format, restricted number, and restricted size.<br />
<br />

# Tech Stack<br />
It is a web application platform developed with the Java,TypeScript and Perl, which is comfortable for users who prefer browser-driven applications or those who are not familiar with UNIX command-line tools.<br />

Task-based user interface design is the most important feature of PlastidHub. The front-end web interface was developed with the vue3 (https://cn.vuejs.org/) framework and Element-Plus (https://element-plus.org/).<br />

The back-end web framework was SpringBoot v3.3.2 (https://docs.spring.io/spring-boot/index.html) with RabbitMQ as message queue and Redis stores result status, which can manage data and parameters submitted by users, and then relay the analysis results (generated files) to the front-end for users to manage (select all, clear selection, download, and delete).<br />

The back-end tools of PlastidHub were developed with Perl v5.26.1 programming language. The dependent third-party tools are BLAST+ v2.16.0, Circos v0.66, MAFFT v7.3, and RAxML v8.2.10.<br />
<br />

# Browser Compatibility<br />
PlastidHub has been tested on recently released web browsers, such as Google Chrome, Microsoft Edge, Mozilla Firefox, Opera, Apple Safari, 360 Speed Browser X, 360 Secure Browser, QQ Browser, and Huawei Browser.<br />
