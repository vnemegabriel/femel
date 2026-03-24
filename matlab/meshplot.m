function meshplot(elementos,nodos,color,showtext)
    patch('Faces',elementos,'Vertices',nodos,'EdgeColor',color,'FaceColor','none', 'LineWidth', 1.5);
    axis equal;
    ax = gca; ax.XTick = []; ax.YTick = []; ax.XColor = 'w'; ax.YColor = 'w';
    if showtext
        for n=1:size(nodos,1)
            text(nodos(n,1),nodos(n,2),[' ' num2str(n)], 'Color', 'blue', 'FontSize', 12, 'FontWeight', 'bold');
        end
        for e=1:size(elementos,1)
            nodes_elem = nodos(elementos(e,:),:);
            text(mean(nodes_elem(:,1)), mean(nodes_elem(:,2)), num2str(e), 'Color','red','FontSize',14,'FontWeight','bold');
        end
    end
end
