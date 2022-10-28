# opt.io - DESIGN

This file documents technical design decisions of opt.io. As **live_optio**, **run_optio**, and all supporting functions are thoroughly commented, many line-by-line considerations are already addressed therein. Here, we will run through some higher-level decisions.

## 1) Modifications to the optimization algorithm

The topology optimization algorithm used herein is based heavily on the work of [Ferrari *et* Sigmund, 2020](https://link.springer.com/article/10.1007/s00158-020-02629-w). The design of the source code, however, offered limited case-by-case flexibility—the material parameters, loads, and constraints were hard-coded, and the only output was the optimization live updates (which MATLAB automatically stitches into a video). Allowing user-specified materials and more convenient load/constraint inputs would provide flexibility, and enabling 3D visualization/mesh exporting would improve the tool's practicality.

Specific, line-by-line implementation of the above is detailed within the **optio.m** comments.

## 2) Separation of live_optio and run_optio

Having two tools for the same job may seem redundant. While I began with **run_optio**, I quickly realized that a beginning user may prefer a few key changes:

1. The ability to quickly change settings without starting over
2. Line-by-line guidance for how to input settings
3. Checking inputs at the final step, as constant error messages can be frustrating
4. A better sense of direction for the inputs and where they lead

This gave birth to **live_optio**. For seasoned users, however, I find that **run_optio** still offers key advantages:

1. A no-frills interface for a more condensed experience
2. Compartmentalized inputs (via independent dialog boxes) to isolate input types and validate inputs frequently (to catch errors as they occur)
3. Easy command-line implementation

## 3) Decision to use MATLAB

Originally, this project was intended to be a web-based app. Time constraints, however, made MATLAB an appealing option to get a working concept off the ground:

1. The original source code is in MATLAB (a 2D Python version exists, but heavy modifications would be required, with math beyond the scope of my understanding).
1. MATLAB is difficult to implement in a web app without dedicated software from MATHWORKS, which I was unable to obtain within the project timeline.
2. With the advanced math required by **optio.m**, the heavy use of MATLAB functions, and differences in indexing, manual transcription of MATLAB to Python was beyond my abilities and available time.
3. With the condensed project timeline, learning how to implement a clean MATLAB GUI would hinder my ability to implement the tool's core functionality. The dialog-box-based and Live-Script-based approaches struck a balance between usability and feasability.

## 4) Misc.

Here, I'll highlight some particularly important, line-level design decisions:

1. Both **live_optio** and **run_optio** automatically clear the workspace before and after they run to keep variables free and the workspace uncluttered.
2. When viewing results, the 3D rotation tool is automatically selected for intuitive 3D orbiting.
3. For **run_optio**, the inability to view all your settings at once means that there is a risk of commiting lots of time & computational power to an misconfigured optimization. The confirmation window before starting the optimization, where users can quickly view their inputs and abort first, handily avoids this issue.
4. The topology optimization can take a long time to run, so accidentally losing one's results can be very frustrating. As such, the "Finish" and "X" (top right of dialog box) buttons at the end of both **live_optio** and **run_optio** will trigger a confirmation message, preventing accidental data losses.
5. The ability to index to the number of elements + 1 may seem odd at first (e.g., 24 x elements allows you to index up to 25 in the x direction). It is, however, an artifact of how the topology optimization is structured (the math of which is beyond my understanding). I felt that preserving this property was important to staying true to the source code and the correct topology optimization algorithm. Artificially capping the user input by removing this +1 "headroom" would work, but it seems to circumvent that which is inherent to the algorithm.