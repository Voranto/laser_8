## Changes
1. Use of a window to average melt depth over the last n iterations. It was called depth_window. Removed it
2. Found a bug in the translation where the logic branching wasn't translated properly. Will be applied in commit "Fixed logic branching bug in translation"
3. In line 454 of the fortran there is some specific branching that checks as follows. If ISTATE > 4 and ISTATE < 8, then it jumps
into a branch that checks if ISTATE == 4. This doesn't make sense as we already ensured that ISTATE > 4. So in my change I moved 
the condition branch back a layer, so the check is independent. This has been removed on ground-truth branch
4. There was a debug ratio I used for a while to test some stuff out. Was currently set to one which equates to
no effect but removed it to make sure.

5. In NVW NVS stuff the translation was I think incorrect, replaced a .GEQ. with a .GT.
6. Some functions were badly ordered? Like the control flow was slightly reordered. Don't really think it had any effect, but 
change it back for consistency sake.

Summary: The translation was quite a mess, the prints should be fixed at some point because they are quite messy but im going to try to first establish a ground truth and then make everything pretty

