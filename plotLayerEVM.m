function plotLayerEVM(NSlots,nslot,pdsch,siz,pdschIndices,pdschSymbols,pdschEq)
% Plot EVM information

persistent slotEVM;
persistent rbEVM
persistent evmPerSlot;

if (nslot==0)
    slotEVM = comm.EVM;
    rbEVM = comm.EVM;
    evmPerSlot = NaN(NSlots,pdsch.NumLayers);
    figure;
end
evmPerSlot(nslot+1,:) = slotEVM(pdschSymbols,pdschEq);
subplot(2,1,1);
plot(0:(NSlots-1),evmPerSlot,'o-');
xlabel('Slot number');
ylabel('EVM (%)');
legend("layer " + (1:pdsch.NumLayers),'Location','EastOutside');
title('EVM per layer per slot');

subplot(2,1,2);
[k,~,p] = ind2sub(siz,pdschIndices);
rbsubs = floor((k-1) / 12);
NRB = siz(1) / 12;
evmPerRB = NaN(NRB,pdsch.NumLayers);
for nu = 1:pdsch.NumLayers
    for rb = unique(rbsubs).'
        this = (rbsubs==rb & p==nu);
        evmPerRB(rb+1,nu) = rbEVM(pdschSymbols(this),pdschEq(this));
    end
end
plot(0:(NRB-1),evmPerRB,'x-');
xlabel('Resource block');
ylabel('EVM (%)');
legend("layer " + (1:pdsch.NumLayers),'Location','EastOutside');
title(['EVM per layer per resource block, slot #' num2str(nslot)]);

drawnow;

end