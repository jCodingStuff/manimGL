<p align="center">
    <a href="https://github.com/3b1b/manim">
        <img src="https://raw.githubusercontent.com/3b1b/manim/master/logo/cropped.png">
    </a>
</p>

[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat)](http://choosealicense.com/licenses/mit/)
[![Manim Subreddit](https://img.shields.io/reddit/subreddit-subscribers/manim.svg?color=ff4301&label=reddit&logo=reddit)](https://www.reddit.com/r/manim/)
[![Manim Discord](https://img.shields.io/discord/581738731934056449.svg?label=discord&logo=discord)](https://discord.com/invite/bYCyhM9Kz2)
[![docs](https://github.com/3b1b/manim/workflows/docs/badge.svg)](https://3b1b.github.io/manim/)

Manim is an engine for precise programmatic animations, designed for creating explanatory math videos.

Note, there are two versions of manim.  This repository began as a personal project by the author of [3Blue1Brown](https://www.3blue1brown.com/) for the purpose of animating those videos, with video-specific code available [here](https://github.com/3b1b/videos).  In 2020 a group of developers forked it into what is now the [community edition](https://github.com/ManimCommunity/manim/), with a goal of being more stable, better tested, quicker to respond to community contributions, and all around friendlier to get started with. See [this page](https://docs.manim.community/en/stable/faq/installation.html#different-versions) for more details.

## Installation
> **WARNING:** These instructions are for ManimGL _only_. Trying to use these instructions to install [ManimCommunity/manim](https://github.com/ManimCommunity/manim) or instructions there to install this version will cause problems. You should first decide which version you wish to install, then only follow the instructions for your desired version.
> 
> **Note**: Unfortunately, you cannot install this fork of `manimgl` through a package manager such as pip. Therefore, you will have to perform a manual installation.

Manim runs on Python 3.7 or higher.

System requirements are [FFmpeg](https://ffmpeg.org/), [OpenGL](https://www.opengl.org/) and [LaTeX](https://www.latex-project.org) (optional, if you want to use LaTeX).
For Linux, [Pango](https://pango.gnome.org) along with its development headers are required. See instruction [here](https://github.com/ManimCommunity/ManimPango#building).

### Directly (Windows)

1. [Install FFmpeg](https://www.wikihow.com/Install-FFmpeg-on-Windows).
2. Install a LaTeX distribution. [MiKTeX](https://miktex.org/download) is recommended.
3. Install the remaining Python packages.
    ```sh
    git clone https://github.com/3b1b/manim.git
    cd manim
    pip install -e .
    manimgl example_scenes.py OpeningManimExample
    ```

### Mac OSX

1. Install FFmpeg, LaTeX in terminal using homebrew.
    ```sh
    brew install ffmpeg mactex
    ```
   
2. Install latest version of manim using these command.
    ```sh
    git clone https://github.com/3b1b/manim.git
    cd manim
    pip install -e .
    manimgl example_scenes.py OpeningManimExample
    ```

## Anaconda Install

1. Install LaTeX as above.
2. Create a conda environment using `conda create -n manim python=3.8`.
3. Activate the environment using `conda activate manim`.
4. Install manimgl using `pip install -e .`.


## Using manim
Try running the following:
```sh
manimgl example_scenes.py OpeningManimExample
```
This should pop up a window playing a simple scene.

Some useful flags include:
* `-w` to write the scene to a file
* `-o` to write the scene to a file and open the result
* `-s` to skip to the end and just show the final frame.
    * `-so` will save the final frame to an image and show it
* `-n <number>` to skip ahead to the `n`'th animation of a scene.
* `-f` to make the playback window fullscreen

Take a look at custom_config.yml for further configuration.  To add your customization, you can either edit this file, or add another file by the same name "custom_config.yml" to whatever directory you are running manim from.  For example [this is the one](https://github.com/3b1b/videos/blob/master/custom_config.yml) for 3blue1brown videos.  There you can specify where videos should be output to, where manim should look for image files and sounds you want to read in, and other defaults regarding style and video quality.

Look through the [example scenes](https://3b1b.github.io/manim/getting_started/example_scenes.html) to get a sense of how it is used, and feel free to look through the code behind [3blue1brown videos](https://github.com/3b1b/videos) for a much larger set of example. Note, however, that developments are often made to the library without considering backwards compatibility with those old videos. To run an old project with a guarantee that it will work, you will have to go back to the commit which completed that project.

### Documentation
Documentation is in progress at [3b1b.github.io/manim](https://3b1b.github.io/manim/). And there is also a Chinese version maintained by [**@manim-kindergarten**](https://manim.org.cn): [docs.manim.org.cn](https://docs.manim.org.cn/) (in Chinese).

[manim-kindergarten](https://github.com/manim-kindergarten/) wrote and collected some useful extra classes and some codes of videos in [manim_sandbox repo](https://github.com/manim-kindergarten/manim_sandbox).


## Contributing
Is always welcome.  As mentioned above, the [community edition](https://github.com/ManimCommunity/manim) has the most active ecosystem for contributions, with testing and continuous integration, but pull requests are welcome here too.  Please explain the motivation for a given change and examples of its effect.


## License
This project falls under the MIT license.


## Contributions of this fork
### Physics engine
Support and animation for a physical system, which is a collection of bodies in space and forces acting among them.

`manimlib.physics.body.Body` represents a physical body with mass, charge, position, and velocity. It can also handle a `manimlib.mobject.mobject.Mobject` instance to represent the body in space and a `manimlib.mobject.geometry.Polyline` instance to track its path through space. The latter is only to be used in 2D animations.

`manimlib.physics.force.Force` is an abstract class represent a force among a set of `Body` instances. It may also manage a set of `Mobject` instances to represent the force graphically. Available subclasses of `Force` are `NewtonGravitationalForce`, `HarmonicBondForce`, `MorseBondForce`, `HarmonicAngleForce`, and `CoulombForce`. Each subclass handles its own parameters. For details, see `manimlib/physics/force.py`.

`manimlib.physics.physical_system.PhysicalSystem` is a subclass of `manimlib.mobject.mobject.Group`, which is in turn a subclass of `Mobject`. This way, `PhysicalSystem` objects can be handled by `manimlib.animation.animation.Animation` instances. It handles a set of `Body` and `Force` instances. There exist subclasses of `PhysicalSystem` such as `GravitationalSystem` that allow the user to just define a set of bodies, and then automatically generate the forces with the `fill_forces` instance method.

`manimlib.animation.physics.EvolvePhysicalSystem` is a subclass of `Animation` that evolves a `PhysicalSystem` instance over time. After the overall force on each body `i` has been computed, it computes the accelerations with `F_i = m_i * a_i`, which are then passed to the time integrator. For details, see `manimlib/animation/physics.py`.<br>
Perhaps an important detail are the `foreground_mobjects` and `background_mobjects` parameters in the `EvolvePhysicalSystem` constructor. Since handling the `Mobject` instances in the system can become complicated, the rendering order of ALL `Mobject` objects in the scene can (and probably will) be altered. Therefore, it is important to let `EvolvePhysicalSystem` know what lays on top and below the `Mobject` instances of the system.

Examples scenes illustrating how to use this engine are located in `example_scenes.py` under the class names `NewtonGravitation2DExample`, `HarmonicBond2DExample`, and `WaterMolecule2DExample`.
