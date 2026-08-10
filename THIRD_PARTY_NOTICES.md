# Third-party notices and clean-room boundary

RDP Web is an independent, MIT-licensed implementation based on published
scientific methods and public user documentation.

No source code from the proprietary Windows RDP package is included. No source
code from PoonLab/OpenRDP is included: OpenRDP is GPL-3.0-licensed. Its bundled
3Seq and GENECONV executables also carry terms that are not suitable for an MIT
distribution, so they are not bundled or invoked here.

AssemblyScript is used only as an Apache-2.0-licensed compiler toolchain. The
compiled WebAssembly module in `public/wasm/rdp.wasm` contains only RDP Web's
clean-room kernels from `assembly/index.ts`.

Scientific method names and citations identify the published procedures that
informed this implementation; they do not imply endorsement by, affiliation
with, or numerical equivalence to the original RDP authors or packages.
