<!--
SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>

SPDX-License-Identifier: GPL-3.0-only
-->

# Contributing to CaΣoS

If you use CaΣoS for your research, [cite us](https://github.com/iFR-OFC/casos#publications) in your publications, report bugs, or extend our code base, you are contributing to CaΣoS and the larger sum-of-squares community! This document details the different ways to contribute to and engage with our repository. To ensure a healthy and respectable community, please adhere to our [Code of Conduct](https://github.com/ifr-ofc/casos#coc-ov-file).

## Questions and Ideas

All questions and ideas should be posted to our [Discussions](https://github.com/iFR-OFC/casos/discussions) boards. Please only create issues if you have a bug or want to suggest a specific enhancement of the software.

## Issues

Please create a [new issue](https://github.com/iFR-OFC/casos/issues/new/choose) if you have encountered a problem with CaΣoS, found a bug, or want to request or suggest a new feature. To help us tackle your issue, be as specific as possible with your description and attach a minimal working example along with full error messages. As far as feature requests are concerned, please be considerate of the fact that CaΣoS is developed and maintained as a research project by academics.

## Code contributions

Any contribution to our codebase must be made via a [pull request](https://github.com/iFR-OFC/casos/compare) to the `main` branch. Pull requests require approval of one of the repository maintainers. Please describe and justify accurately all your changes and additions. All code must adhere to our coding guidelines and be compatible with our licence. Major changes, including new features and interfaces, should be coordinated in advance to ensure that we provide adequate support.

### Coding guidelines

All code contributions must follow the [CaΣoS Coding Guidelines](https://github.com/iFR-OFC/casos-coding-guidelines) which are based on the MATLAB Coding Guidelines.

### Copyright & License

It is understood that the copyright for all code is given to the CaΣoS maintainers. We use the [REUSE specification](https://reuse.software/spec-3.3/) for copyright declaration and organisation of license declarations in all files included in the repository. All new files therefore must contain the following header in a comment of the respective programming languages:

<!-- REUSE-IgnoreStart -->
```matlab
% SPDX-FileCopyrightText: YEAR Institute of Flight Mechanics and Controls, Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and FURTHER AUTHORS <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileCopyrightText: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only
```
<!-- REUSE-IgnoreEnd -->

Contributors may become authors of the respective files if they are responsible for a significant part of the addition. Other changes to either copyright or the [license](https://github.com/iFR-OFC/casos#GPL-3.0-1-ov-file) are generally not permitted. In exceptional cases, third-party code can be included under a different license contingent on compliance. In this case, the relevant inclusions including details of the copyright holder(s), version history, and license must be identified using REUSE's [in-line snippet comments](https://reuse.software/spec-3.3/#in-line-snippet-comments) as follows:

<!-- REUSE-IgnoreStart -->
```matlab
% SPDX-SnippetBegin
% SPDX-SnippetCopyrightText: 2022 Jane Doe <jane@example.com>
% SPDX-License-Identifier: MIT

print("Hello, world!");

% SPDX-SnippetEnd
```
<!-- REUSE-IgnoreEnd -->

## Examples

You are always welcome to add examples to the [CaΣoS example package](https://github.com/iFR-OFC/casos-example-package). We are particularly interested in innovative ways to use sum-of-squares optimization and their applications to systems and control problems from various disciplines.

## Academic cooperations

We are always open to academic cooperations and involvement in research projects. All inquiries should come with a clear project proposal and publication perspective.

----

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://github.com/ifr-ofc/.github/blob/main/assets/logo-casos-inverted.png">
  <source media="(prefers-color-scheme: light)" srcset="https://github.com/ifr-ofc/.github/blob/main/assets/logo-casos-trans.png">
  <img width="50%" alt="CaΣoS: A nonlinear sum-of-squares optimization suite">
</picture>
