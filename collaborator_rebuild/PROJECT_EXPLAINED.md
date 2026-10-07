# The SURD project in simple terms

We are trying to understand how the changing light from an active galaxy
carries information about the gas around its central black hole. We want to
learn which measurements tell us something new, which repeat the same
information, and which become useful only when considered together.

The galaxy is NGC 5548. Near its central black hole, hot material produces
bright, variable radiation. Surrounding gas absorbs some of that radiation
and emits light of its own. One emission feature is called H beta, produced
by hydrogen.

When the central source brightens, the gas can respond later because light
takes time to travel. Studying these delayed responses is called reverberation
mapping. It helps astronomers investigate the size and motion of gas near the
black hole. [Peterson and Wandel](https://arxiv.org/abs/astro-ph/0007147)

Think of a lamp surrounded by material that glows when illuminated. The lamp
changes brightness, and the surrounding material responds. Our observations
are a distant, unresolved record of those changes.

We use five brightness records. Each record is a list of dates and measured
brightnesses, often called a light curve.

| Measurement | Simple meaning |
| :--- | :--- |
| F5100 continuum | Optical brightness near the central source |
| Total H beta | Brightness of the whole broad hydrogen emission feature |
| Blue H beta | Brightness in its shorter wavelength side |
| Core H beta | Brightness in its central part |
| Red H beta | Brightness in its longer wavelength side |

The different parts of the emission feature provide information about
different projected gas velocities. They are not necessarily three separate
physical clouds. Also, the optical continuum is a proxy for the radiation
illuminating the gas.
[Cackett and Horne](https://academic.oup.com/mnras/article/365/4/1180/992422)

Our specific question is about information. For example:

> How much does knowing the continuum and the different line components ten
> days earlier help us predict the blue component today?

We also include the blue component's own earlier brightness. Otherwise, its
natural persistence could be mistaken for information supplied by another
measurement. We repeat the calculation with different targets and delays.

The theory we use is SURD, which divides predictive information into three
kinds:

| Kind | Simple meaning | Illustrative example |
| :--- | :--- | :--- |
| Unique | Information available from one measurement alone | Blue contains a clue the other measurements do not provide |
| Redundant | Information shared by several measurements | Blue and red both repeat the same clue |
| Synergistic | Information available through a combination | Blue and red together reveal something neither reveals alone |

For synergy, imagine two random switches. A light turns on when their
positions differ. Knowing either switch alone does not tell you whether the
light is on. Knowing both does. That is an example of information carried by
a combination.

SURD also measures leakage, meaning the uncertainty about the target that
remains after using the predictors. Its foundation is information theory:
useful information reduces uncertainty about what happens next. Applying
this to astronomical observations requires careful checks before interpreting
a predictive relationship as a physical cause.
[The SURD paper](https://www.nature.com/articles/s41467-024-53373-4)

What we are doing now is checking whether our measurements can support that
analysis reliably.

The calculation estimates probabilities by counting how often different
brightness combinations occur. For our current method, we divide brightness
into categories:

1. With two bins, the categories are lower and higher brightness.
2. With three bins, they are lower, middle, and higher brightness.

With one target and four predictors, there are 32 possible category
combinations using two bins, and 243 using three. We often have only tens of
usable observations. Many combinations therefore occur rarely or never,
making their estimated probabilities uncertain.

We have been checking two things. First, does the answer change substantially
when we use two categories instead of three? Second, does it change because
different delays use different observation dates?

To address the second question, we recently required every delay from 1
through 30 days to use exactly the same target dates within each curve.
That makes the comparison fairer, but removes many observations.

Only four curves retained enough observations to pass our existing numerical
guard. All four predict the continuum from its own history and past line
components. Their information breakdown still changes substantially with bin
choice. The current evidence therefore does not establish a stable physical
delay. The details are in our [fixed support review](FIXED_SUPPORT_REVIEW.md).
This status describes the calculations completed on 7 October 2026.

The vision is to build a trustworthy map of how information moves through
these observed signals. We hope to identify whether the line components
mainly share a common brightness variation, whether some carry additional
information, and whether combinations reveal relationships that individual
measurements miss.

Our next step is calibration: test the calculation on artificial data where
we already know the correct relationships, using observation counts comparable
to ours. That will help establish what the method can recover reliably.
A scientifically useful outcome could be either a robust astronomical
relationship or a clear demonstration that these historical observations
cannot support this particular decomposition.
