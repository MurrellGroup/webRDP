# Third-party notices and source provenance

RDP Web is an MIT-licensed browser implementation based on published methods,
the RDP5 manual, and author-supplied RDP5 source.

The original RDP authors supplied the RDP5 Visual Basic/native source and gave
permission to use it for this project. Source-compatible ports of selected
active routines are therefore included, with the legacy routine names recorded
in `SOURCE_WORKFLOW_AUDIT.md` and `SOURCE_PARITY.md`. The VB6/Win32 application,
DLL binaries, and source archives themselves are not redistributed.

The RDP5 source identifies its `PHITest2` helper family as a Visual Basic
translation of Trevor Bruen's PHIPACK implementation. RDP Web ports the
observable RDP5 formulas and graph algorithm from the supplied source and
attributes the PHI method to Bruen, Philippe & Bryant (2006); it does not bundle
PHIPACK source or binaries.

No source code from PoonLab/OpenRDP is included; OpenRDP is GPL-3.0-licensed.
Its bundled external executables are also not bundled or invoked here.

AssemblyScript is used only as an Apache-2.0-licensed compiler toolchain. The
compiled WebAssembly module in `public/wasm/rdp.wasm` contains RDP Web kernels,
including the permitted RDP5 source-compatibility ports in `assembly/index.ts`.

Scientific method names and citations identify the procedures being
implemented. A port is labeled source-compatible only where its active source
path has been mapped and tested; that label is not a blanket claim of complete
desktop numerical parity.
