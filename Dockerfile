FROM --platform=linux/amd64 mcr.microsoft.com/dotnet/runtime:8.0
COPY MaxQuant_v.2.9.0 /opt/maxquant
WORKDIR /work
ENTRYPOINT ["dotnet", "/opt/maxquant/bin/MaxQuantCmd.dll"]
