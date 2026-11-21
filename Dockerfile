FROM --platform=linux/amd64 mcr.microsoft.com/dotnet/runtime:8.0
# libgomp1 (OpenMP): needed by MaxQuant's libxgboost.so (MaxDIA library prediction) and
# libbaf2sql_c.so (Bruker raw files); the .NET base image lacks it
RUN apt-get update \
    && apt-get install -y --no-install-recommends libgomp1 \
    && rm -rf /var/lib/apt/lists/*
COPY MaxQuant_v.2.9.0 /opt/maxquant
WORKDIR /work
ENTRYPOINT ["dotnet", "/opt/maxquant/bin/MaxQuantCmd.dll"]
